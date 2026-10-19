"""Scheduler-free online output with real local MPI and the production writer.

Run: python3 tests/test_output_ordered_online.py -v
MPICC/MPIEXEC may select a local MPI installation. All artifacts are temporary.
"""
import ctypes
import ctypes.util
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest
import zlib

from test_output_discard import MASK, reference_hash
from test_sam_stream_md5 import build_md5_object

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ((0, 1), (0, 6), (1, 6))


def blobs(pass_id, chunk, zero_bytes=False):
    reads = 137 - pass_id * 30 if chunk % 3 == 0 else 5 + chunk % 5 - pass_id
    result = []
    for read in range(reads):
        lines = []
        if not (zero_bytes or chunk % 4 == 1 or read % 7 == 6):
            for record in range(1 + (chunk + read) % 3):
                line = (f"p{pass_id}-c{chunk:03d}-r{read:03d}\t{2048 if record else 0}"
                        f"\tchr{record + 1}\t{chunk * 1000 + read + 1}"
                        "\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII").encode()
                if read == 100 and record == 0:
                    line += b"\tZZ:Z:" + bytes([ord("A") + chunk % 20]) * (70000 - pass_id * 123)
                    line += bytes([ord("a") + pass_id])
                lines.append(line + b"\n")
        result.append(b"".join(lines))
    return result


def reports(stderr, label):
    pattern = rf"\[SWBWA {re.escape(label)} rank (\d+)/(\d+)\] ([^\n]*)"
    return [(int(rank), int(size),
             {key: int(value, 16) if value.startswith("0x") else int(value)
              for key, value in (item.split("=") for item in line.split())})
            for rank, size, line in re.findall(pattern, stderr)]


class OrderedOnlineOutputTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="swbwa-ordered-online-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory = Path(cls.temporary.name)
        cls.executables = {}
        library = ctypes.util.find_library("z")
        if library is None:
            raise RuntimeError("libz is required for independent ordered CRC combination")
        cls.libz = ctypes.CDLL(library)
        cls.crc_combine = cls.libz.crc32_combine
        cls.crc_combine.argtypes = (ctypes.c_ulong, ctypes.c_ulong, ctypes.c_long)
        cls.crc_combine.restype = ctypes.c_ulong
        md5_link = build_md5_object(cls.directory)
        for only, workers in CONFIGS:
            executable = cls.directory / f"online-{only}-{workers}"
            subprocess.run(shlex.split(os.getenv("MPICC", "mpicc")) + [
                "-std=gnu11", "-O2", "-Wall", "-Wextra", "-Werror", "-Wno-unused-function",
                "-Iinclude", "-DSWBWA_USE_MPI=1", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                "-DSWBWA_ENABLE_GZIP_INPUT=1", "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_STATIC",
                "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_ORDERED",
                "-DSWBWA_OUTPUT_BUFFER_BYTES=16", f"-DSWBWA_OUTPUT_RMA_ONLY={only}",
                f"-DSWBWA_HOST_MPE_THREADS={workers}", "tests/test_output_ordered_online.c",
                "src/host/swbwa_mpi.c", "src/host/swbwa_output.c",
                "src/host/swbwa_host_workers.c", "-pthread", "-o", str(executable),
            ] + md5_link, cwd=ROOT, check=True, timeout=60)
            cls.executables[only, workers] = executable

    def run_case(self, only, workers, ranks, chunks, mode="roundrobin", sample_limit=0):
        path = self.directory / f"out-{only}-{workers}-{ranks}-{chunks}-{mode}-{sample_limit}.sam"
        env = dict(os.environ)
        env.pop("SWBWA_OUTPUT_MD5", None)
        env.pop("SWBWA_OUTPUT_SAMPLE_READS", None)
        # Keep Hydra on the local machine even in a cluster login environment.
        env["HYDRA_LAUNCHER"] = "fork"
        env.pop("HYDRA_HOST_FILE", None)
        if sample_limit is not None:
            env["SWBWA_OUTPUT_SAMPLE_READS"] = str(sample_limit)
        command = shlex.split(os.getenv("MPIEXEC", "mpiexec")) + [
            "-n", str(ranks), str(self.executables[only, workers]), str(path), str(chunks), mode,
        ]
        try:
            run = subprocess.run(command, capture_output=True, text=True, env=env,
                                 cwd=ROOT, timeout=60)
        except subprocess.TimeoutExpired as error:
            self.fail(f"local MPI timed out: {command!r}\n{error.stdout!r}\n{error.stderr!r}")
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertIn("online ordered output PASS", run.stdout)
        zero_bytes = mode == "zero-bytes"
        expected = b"".join(blob for chunk in range(chunks) for blob in blobs(1, chunk, zero_bytes))
        if only:
            self.assertFalse(path.exists(), "RMA-only mode created an output file")
        else:
            self.assertEqual(path.read_bytes(), expected)
        for pass_id in range(2):
            blocks = re.findall(rf"\[online test\] pass={pass_id} begin\n(.*?)"
                                rf"\[online test\] pass={pass_id} end\n", run.stderr, re.S)
            self.assertEqual(len(blocks), 1, run.stderr)
            self.check_reports(blocks[0], pass_id, only, ranks, chunks, zero_bytes,
                               100 if sample_limit is None else sample_limit)
        return run

    def check_reports(self, stderr, pass_id, only, ranks, chunks, zero_bytes, sample_limit):
        expected = {chunk: blobs(pass_id, chunk, zero_bytes) for chunk in range(chunks)}
        extents = reports(stderr, "output extent")
        self.assertEqual(len(extents), chunks, stderr)
        self.assertEqual(sorted(entry[2]["chunk"] for entry in extents), list(range(chunks)))
        offset = 0
        for rank, size, entry in sorted(extents, key=lambda item: item[2]["chunk"]):
            chunk = entry["chunk"]
            self.assertEqual((rank, size), (chunk % ranks, ranks))
            payload = b"".join(expected[chunk])
            self.assertEqual(entry["start"], chunk * 10000)
            self.assertEqual(entry["end"], (chunk + 1) * 10000)
            self.assertEqual(entry["reads"], len(expected[chunk]))
            self.assertEqual(entry["bytes"], len(payload))
            self.assertEqual(entry["offset"], offset)
            offset += len(payload)
        self.assertEqual(offset, sum(len(blob) for values in expected.values() for blob in values))
        samples = reports(stderr, "chunk sample")
        summaries = reports(stderr, "output RMA-only")
        if not only:
            self.assertEqual(samples, [])
            self.assertEqual(summaries, [])
            return
        self.assertEqual(len(samples), chunks, stderr)
        self.assertEqual(sorted(entry[2]["chunk"] for entry in samples), list(range(chunks)))
        for rank, size, sample in samples:
            chunk = sample["chunk"]
            self.assertEqual((rank, size), (chunk % ranks, ranks))
            values = expected[chunk]
            sampled = values if sample_limit == 0 else values[:sample_limit]
            hashes = [reference_hash(blob) for blob in sampled]
            xor = 0
            for value in hashes:
                xor ^= value
            self.assertEqual(sample["start"], chunk * 10000)
            self.assertEqual(sample["end"], (chunk + 1) * 10000)
            self.assertEqual(sample["paired"], 0)
            self.assertEqual(sample["reads"], len(values))
            self.assertEqual(sample["sam_bytes"], sum(map(len, values)))
            self.assertEqual(sample["sample_reads"], len(sampled))
            self.assertEqual(sample["sample_bytes"], sum(map(len, sampled)))
            self.assertEqual(sample["sum"], sum(hashes) & MASK)
            self.assertEqual(sample["xor"], xor)
            self.assertEqual(sample["sam_records"], sum(blob.count(b"\n") for blob in sampled))
            if sample_limit == 0:
                self.assertEqual(sample["packed_crc32"], zlib.crc32(b"".join(values)))
        if sample_limit == 0:
            packed = {sample["chunk"]: sample for _, _, sample in samples}
            combined = 0
            for _, _, extent in sorted(extents, key=lambda item: (item[2]["offset"], item[2]["chunk"])):
                sample = packed[extent["chunk"]]
                combined = self.crc_combine(combined, sample["packed_crc32"], extent["bytes"])
            ordered = b"".join(blob for chunk in range(chunks) for blob in expected[chunk])
            self.assertEqual(combined, zlib.crc32(ordered))
        self.assertEqual(len(summaries), ranks, stderr)
        self.assertEqual(sorted(entry[0] for entry in summaries), list(range(ranks)))
        for rank, size, summary in summaries:
            self.assertEqual(size, ranks)
            owned = list(range(rank, chunks, ranks))
            values = [blob for chunk in owned for blob in expected[chunk]]
            self.assertEqual(summary["chunks"], len(owned))
            self.assertEqual(summary["reads"], len(values))
            self.assertEqual(summary["sam_bytes"], sum(map(len, values)))
            self.assertEqual(summary["reserved_bytes"], sum(map(len, values)))
            self.assertEqual(summary["reservations"], len(owned))
            self.assertEqual(summary["disk_bytes"], 0)
            self.assertEqual(summary["sample_limit"], sample_limit)

    def test_ordered_bytes_full_hash_and_reopen(self):
        for only, workers in CONFIGS:
            for ranks in (1, 3, 4):
                cases = [(0, "roundrobin"), (1, "roundrobin"),
                         (2 * ranks + 3, "roundrobin"), (ranks + 2, "zero-bytes"),
                         (2 * ranks + 3, "decreasing")]
                if ranks > 1:
                    cases.append((ranks - 1, "roundrobin"))
                for chunks, mode in cases:
                    with self.subTest(only=only, workers=workers, ranks=ranks, chunks=chunks, mode=mode):
                        run = self.run_case(only, workers, ranks, chunks, mode)
                        self.assertEqual(run.stdout.count("online output bytes pass="), 2)
                        if ranks > 1 and chunks < ranks:
                            self.assertEqual(run.stdout.count("online local finish PASS"), 2)
                        if mode == "decreasing":
                            self.assertEqual(run.stdout.count("online decreasing ID rejected PASS"), 2)

    def test_close_rejects_missing_or_wrong_finish(self):
        for only, workers in CONFIGS:
            for ranks in (1, 3, 4):
                for mode, chunks in (("missing-finish", 0), ("missing-finish", ranks + 2),
                                     ("missing-one-finish", 1), ("wrong-low", ranks + 2),
                                     ("wrong-high", 0), ("inconsistent-count", ranks + 2)):
                    with self.subTest(only=only, workers=workers, ranks=ranks, chunks=chunks, mode=mode):
                        run = self.run_case(only, workers, ranks, chunks, mode)
                        self.assertEqual(run.stdout.count("online close error pass="), 2)

    def test_default_sampling_remains_100(self):
        for ranks in (1, 3, 4):
            with self.subTest(ranks=ranks):
                self.run_case(1, 6, ranks, ranks + 2, sample_limit=None)


if __name__ == "__main__":
    unittest.main()
