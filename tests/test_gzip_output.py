#!/usr/bin/env python3
"""Exercise production scatter/gather gzip output with native zlib and MPI."""
import gzip
from functools import lru_cache
import os
import re
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
import zlib

from test_sam_stream_md5 import build_md5_object

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(maxsize=1)
def random_ascii_body(length):
    body = bytearray(length)
    for i in range(length):
        mixed = (i * 2654435761) & 0xffffffff
        mixed ^= mixed >> 13
        mixed = (mixed * 0x85ebca6b) & 0xffffffff
        mixed ^= mixed >> 16
        body[i] = 33 + mixed % 90
    return bytes(body)


def expected(chunks):
    values = []
    alphabet = b"ACGTN0123456789"
    for chunk in chunks:
        for read in range(31):
            length = (read * 73 + chunk * 41) % 4000
            if read == 5 and chunk % 2 == 0:
                length += 1048576
            if read % 7 == 0:
                continue
            payload = bytes(alphabet[(i * 7 + read + chunk) % 15] for i in range(15))
            body = (payload * (length // 15 + 1))[:length]
            if chunk == 2 and read == 5:
                body = random_ascii_body(length)
            values.append(f"chunk{chunk:03d}-read{read:03d}-pass1\t".encode() + body + b"\n")
    return b"".join(values)


class GzipOutputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="swbwa-gzip-output-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.binaries = {}
        for mpi, mode, workers in ((0, "split", 1), (0, "split", 6),
                                    (1, "split", 1), (1, "single_unordered", 1),
                                    (1, "single_ordered", 6)):
            name = f"{mpi}-{mode}-{workers}"
            binary = cls.directory / name
            compiler = shlex.split(os.getenv("MPICC", "mpicc"))
            command = compiler + ["-std=gnu11", "-O2", "-g", "-Wall", "-Wextra",
                "-Werror", "-Wno-unused-function", "-D_GNU_SOURCE", "-D_FILE_OFFSET_BITS=64",
                "-Iinclude", "-Isrc/gzip", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                "-DSWBWA_ENABLE_GZIP_OUTPUT=1", "-DSWBWA_ENABLE_GZIP_INPUT=1",
                f"-DSWBWA_USE_MPI={mpi}", f"-DSWBWA_HOST_MPE_THREADS={workers}",
                "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_DYNAMIC",
                f"-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_{mode.upper()}",
                "tests/test_gzip_output.c", "src/host/swbwa_output.c", "src/host/swbwa_mpi.c",
                "src/host/swbwa_host_workers.c", "src/gzip/output.c", "src/gzip/codec_job.c",
                "src/gzip/backend_host.c", "-pthread", "-lz", "-o", str(binary)]
            command += ["tests/gzip_output_short_write.c",
                        "-Wl,--wrap=writev,--wrap=pwritev64,--wrap=swbwa_codec_launch,--wrap=sysconf",
                        "-Wl,--wrap=swbwa_mpi_fastq_scheduler_chunk_count,--wrap=swbwa_mpi_fastq_scheduler_ticket_mode"]
            command += build_md5_object(cls.directory)
            sanitizer = os.environ.get("SWBWA_TEST_SANITIZER")
            if sanitizer:
                command += ["-fsanitize=" + sanitizer, "-fno-sanitize-recover=all",
                            "-fno-omit-frame-pointer", "-fno-pie", "-no-pie"]
            subprocess.run(command, cwd=ROOT, check=True)
            cls.binaries[mpi, mode, workers] = binary
            if mpi and mode == "single_ordered":
                cls.rma_binary = cls.directory / "ordered-rma-only"
                command[command.index(str(binary))] = str(cls.rma_binary)
                command.append("-DSWBWA_OUTPUT_RMA_ONLY=1")
                subprocess.run(command, cwd=ROOT, check=True)

    def run_case(self, mpi, mode, workers, ranks, chunks, block=0, level=1, suffix=".sam.gz", short=False, no_pool=False, parallel=False,
                 prepare=False, api=False, slow=False, empty=False, errors=False, poll_us=None,
                 prefix_ahead=False, known_prefix=False, iov_limit=None):
        path = self.directory / f"out-{mpi}-{mode}-{workers}-{ranks}-{chunks}-{block}-{level}{suffix}"
        command = [str(self.binaries[mpi, mode, workers]), str(path), str(chunks)]
        if mpi:
            command = shlex.split(os.getenv("MPIEXEC", "mpiexec")) + ["-n", str(ranks)] + command
        env = dict(os.environ, HYDRA_LAUNCHER="fork", SWBWA_GZIP_BLOCK_KIB=str(block),
                   SWBWA_GZIP_LEVEL=str(level))
        env.pop("SWBWA_OUTPUT_MD5", None)
        env.pop("HYDRA_HOST_FILE", None)
        env.pop("SWBWA_TEST_IOV_LIMIT", None)
        if iov_limit is not None:
            env["SWBWA_TEST_IOV_LIMIT"] = str(iov_limit)
        env["SWBWA_GZIP_PREFIX_AHEAD"] = str(int(prefix_ahead))
        env.pop("SWBWA_TEST_KNOWN_PREFIX", None)
        if known_prefix:
            env["SWBWA_TEST_KNOWN_PREFIX"] = "1"
        if poll_us is not None:
            env["SWBWA_ORDERED_POLL_US"] = str(poll_us)
        for key in ("SWBWA_TEST_PREPARE_API", "SWBWA_TEST_SLOW_WRITER", "SWBWA_TEST_EMPTY_RECORDS",
                    "SWBWA_TEST_POOL_ERRORS", "SWBWA_TEST_WRITE_FAILURE", "SWBWA_TEST_NO_HOST_POOL",
                    "SWBWA_TEST_WRITE_LIMIT", "SWBWA_TEST_WRITE_EINTR"):
            env.pop(key, None)
        if prepare is None:
            env.pop("SWBWA_GZIP_PREPARE_OUTPUT", None)
        else:
            env["SWBWA_GZIP_PREPARE_OUTPUT"] = "1" if prepare else "0"
        if api or prepare:
            env["SWBWA_TEST_PREPARE_API"] = "1"
        if slow:
            env["SWBWA_TEST_SLOW_WRITER"] = "1"
        if empty:
            env["SWBWA_TEST_EMPTY_RECORDS"] = "1"
        if errors:
            env["SWBWA_TEST_POOL_ERRORS"] = "1"
        env.setdefault("ASAN_OPTIONS", "detect_leaks=0")
        if short:
            env.update(SWBWA_TEST_WRITE_LIMIT="777", SWBWA_TEST_WRITE_EINTR="1")
        if no_pool:
            env["SWBWA_TEST_NO_HOST_POOL"] = "1"
        env["SWBWA_GZIP_PARALLEL_WRITE"] = "1" if parallel else "0"
        run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn("gzip output driver PASS", run.stdout)
        outputs = [(path, range(chunks))]
        if mpi and mode == "split":
            outputs = [(path.with_name(path.name[:-len(suffix)] + f".rank{rank:06d}" + suffix),
                        range(rank, chunks, ranks)) for rank in range(ranks)]
        for output, owned in outputs:
            encoded = output.read_bytes()
            decoded = gzip.decompress(encoded) if suffix.endswith(".gz") else encoded
            reference = b"" if empty else expected(owned)
            if mode == "single_unordered":
                self.assertEqual(sorted(decoded.splitlines()), sorted(reference.splitlines()))
            else:
                self.assertEqual(decoded, reference)
            if suffix.endswith(".gz"):
                subprocess.run(["gzip", "-t", str(output)], check=True)
                members = 0
                while encoded:
                    stream = zlib.decompressobj(31)
                    stream.decompress(encoded)
                    self.assertTrue(stream.eof)
                    encoded = stream.unused_data
                    members += 1
                self.assertGreater(members, 0)
        if suffix.endswith(".gz") and not errors:
            if iov_limit is not None:
                self.assertIn(f"iov_limit={iov_limit if iov_limit > 0 else 16}", run.stderr)
            self.assertIn("gzip write wall bandwidth", run.stderr)
            self.assertIn("Parallel syscall times are worker sums, not wall time.", run.stderr)
            reports = re.findall(r"\[SWBWA gzip output rank [^]]+\][^\n]*batches=(\d+)[^\n]*raw_bytes=(\d+)", run.stderr)
            self.assertEqual(len(reports), 2 * ranks)
            self.assertEqual(sum(int(batches) for batches, _ in reports), 2 * chunks)
            self.assertEqual(sum(int(raw) for _, raw in reports), 0 if empty else 2 * len(expected(range(chunks))))
        return run

    def test_modes(self):
        for mpi, mode, workers in self.binaries:
            for ranks, chunks in ((1, 7), (3, 7), (3, 1), (3, 0)) if mpi else ((1, 7), (1, 0)):
                with self.subTest(mpi=mpi, mode=mode, workers=workers, ranks=ranks, chunks=chunks):
                    self.run_case(mpi, mode, workers, ranks, chunks)

    def test_parameters(self):
        for block, level in ((1, 1), (64, 4), (4096, 9)):
            with self.subTest(block=block, level=level):
                self.run_case(0, "split", 6, 1, 3, block, level)

    def test_short_writes(self):
        self.run_case(0, "split", 1, 1, 3, block=1, short=True)
        self.run_case(1, "single_ordered", 6, 3, 7, block=1, short=True, parallel=True)

    def test_iov_limits(self):
        for limit in (-1, 16, 128, 1024):
            with self.subTest(limit=limit):
                self.run_case(1, "single_ordered", 6, 3, 7, block=1,
                              prepare=True, prefix_ahead=True, iov_limit=limit)

    def test_pool_not_initialized(self):
        self.run_case(0, "split", 6, 1, 3, no_pool=True, parallel=True)

    def test_parallel_writes(self):
        self.run_case(0, "split", 6, 1, 3, parallel=True)
        self.run_case(1, "single_ordered", 6, 3, 7, parallel=True)

    def test_ordered_poll_pause(self):
        for pause in (0, 1000):
            with self.subTest(poll_us=pause):
                self.run_case(1, "single_ordered", 6, 3, 7, parallel=True,
                              prepare=True, slow=True, poll_us=pause)

    def test_prefix_ahead(self):
        for known, chunks, empty in ((False, 19, False), (True, 19, False),
                                     (False, 7, True), (True, 0, False)):
            with self.subTest(known_prefix=known, chunks=chunks, empty=empty):
                result = self.run_case(1, "single_ordered", 6, 3, chunks, prepare=True,
                                      slow=True, parallel=True, empty=empty,
                                      prefix_ahead=True, known_prefix=known)
                self.assertEqual(len(re.findall(r'\[SWBWA gzip prefix rank ', result.stderr)), 6)

    def test_prepare_modes_and_slow_writer(self):
        for mpi, mode, workers in self.binaries:
            with self.subTest(mpi=mpi, mode=mode, workers=workers):
                self.run_case(mpi, mode, workers, 3 if mpi else 1, 19,
                              prepare=True, slow=True, parallel=True)

    def test_default_prepares_on_producer_with_slow_writer(self):
        for mpi, mode, workers in self.binaries:
            with self.subTest(mpi=mpi, mode=mode, workers=workers):
                self.run_case(mpi, mode, workers, 3 if mpi else 1, 19,
                              prepare=None, api=True, slow=True, parallel=True)

    def test_prepare_empty(self):
        self.run_case(0, "split", 6, 1, 7, prepare=True, slow=True, empty=True)
        self.run_case(1, "single_ordered", 6, 3, 7, prepare=True, slow=True, empty=True)
        self.run_case(1, "single_ordered", 6, 3, 0, prepare=True, slow=True)

    def test_prepare_fallbacks(self):
        self.run_case(0, "split", 1, 1, 7, api=True, prepare=False)
        self.run_case(0, "split", 6, 1, 7, suffix=".sam", api=True, prepare=None, slow=True)
        self.run_case(1, "single_ordered", 6, 3, 7, suffix=".sam", api=True, prepare=None, slow=True)
        self.run_case(0, "split", 6, 1, 7, suffix=".sam", prepare=True, slow=True)
        self.run_case(1, "single_ordered", 6, 3, 7, suffix=".sam", prepare=True, slow=True)
        self.run_case(1, "split", 1, 3, 7, suffix=".sam", prepare=True, slow=True)

    def test_prepare_short_writes(self):
        self.run_case(0, "split", 1, 1, 7, block=1, short=True, prepare=True, slow=True)
        self.run_case(1, "single_ordered", 6, 3, 19, block=1, short=True,
                      prepare=True, slow=True, parallel=True)

    def test_prepare_pool_errors(self):
        self.run_case(0, "split", 6, 1, 7, prepare=True, errors=True)

    def test_prepare_without_host_pool(self):
        self.run_case(0, "split", 6, 1, 7, prepare=True, slow=True, no_pool=True, parallel=True)

    def test_prepare_write_error_releases_slots(self):
        path = self.directory / "full.sam.gz"
        path.symlink_to("/dev/full")
        env = dict(os.environ, SWBWA_GZIP_PREPARE_OUTPUT="1", SWBWA_TEST_PREPARE_API="1",
                   SWBWA_TEST_WRITE_FAILURE="1", SWBWA_GZIP_PARALLEL_WRITE="0")
        env.pop("SWBWA_OUTPUT_MD5", None)
        env.setdefault("ASAN_OPTIONS", "detect_leaks=0")
        run = subprocess.run([str(self.binaries[0, "split", 1]), str(path), "0"],
                             cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn("gzip output driver PASS", run.stdout)

    def run_rma_only(self, prepare):
        path = self.directory / "rma-only.sam.gz"
        env = dict(os.environ, HYDRA_LAUNCHER="fork", SWBWA_OUTPUT_SAMPLE_READS="0")
        if prepare is None:
            env.pop("SWBWA_GZIP_PREPARE_OUTPUT", None)
        else:
            env["SWBWA_GZIP_PREPARE_OUTPUT"] = "1" if prepare else "0"
        env["SWBWA_TEST_PREPARE_API"] = "1"
        if prepare is not False:
            env["SWBWA_TEST_SLOW_WRITER"] = "1"
        env.setdefault("ASAN_OPTIONS", "detect_leaks=0")
        command = shlex.split(os.getenv("MPIEXEC", "mpiexec")) + [
            "-n", "3", str(self.rma_binary), str(path), "7"]
        run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertFalse(path.exists())
        checksums = dict((int(chunk), int(crc, 16)) for chunk, crc in re.findall(
            r"\[SWBWA chunk sample rank [^]]+\] chunk=(\d+)[^\n]*packed_crc32=0x([a-f0-9]+)",
            run.stderr))
        self.assertEqual(set(checksums), set(range(7)), run.stderr)
        for chunk, crc in checksums.items():
            self.assertEqual(crc, zlib.crc32(expected(range(chunk, chunk + 1))))

    def test_rma_only_compresses_and_hashes_original_sam(self):
        self.run_rma_only(False)

    def test_prepared_rma_only(self):
        self.run_rma_only(True)

    def test_default_prepared_rma_only(self):
        self.run_rma_only(None)


if __name__ == "__main__":
    unittest.main()
