"""Real-MPI writer test; standalone, needs mpicc/mpiexec and three local ranks."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

from test_output_discard import reference_hash, MASK
from test_sam_stream_md5 import build_md5_object

ROOT = Path(__file__).resolve().parents[1]


def records(rank, chunk):
    result = []
    for i in range(2 if chunk else 104):
        data = f"r{rank}-c{chunk}-i{i}\t".encode()
        data += bytes([ord("A") + i % 20]) * (257 if i == 0 else 1 + i % 23)
        data += b"\n"
        if i == 0:
            data += b"supplementary\n"
        result.append(data)
    return result


class OutputRmaTest(unittest.TestCase):
    def test_chunk_extents_and_parallel_packing(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            md5_link = build_md5_object(directory)
            for threads in (1, 6):
                for only in (0, 1):
                    with self.subTest(threads=threads, rma_only=only):
                        exe = directory / f"chunks-{threads}-{only}"
                        subprocess.run(shlex.split(os.environ.get("MPICC", "mpicc")) + [
                            "-std=gnu11", "-Wall", "-Wextra", "-Werror", "-Wno-unused-function",
                            "-Iinclude", "-DSWBWA_USE_MPI=1", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                            "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_UNORDERED",
                            "-DSWBWA_OUTPUT_BUFFER_BYTES=128", f"-DSWBWA_OUTPUT_RMA_ONLY={only}",
                            f"-DSWBWA_HOST_MPE_THREADS={threads}", "-DTEST_CHUNK_OUTPUT=1",
                            "tests/test_output_rma.c", "src/host/swbwa_output.c",
                            "src/host/swbwa_host_workers.c", "-pthread", "-o", str(exe)] + md5_link,
                            cwd=ROOT, check=True)
                        path = directory / f"chunks-{threads}-{only}.sam"
                        run = subprocess.run(shlex.split(os.environ.get("MPIEXEC", "mpiexec")) + [
                            "-n", "3", str(exe), str(path)], capture_output=True, text=True)
                        self.assertEqual(run.returncode, 0, run.stderr)
                        entries = re.findall(r"\[SWBWA output extent rank (\d+)/(\d+)\] (.*)", run.stderr)
                        self.assertEqual(len(entries), 4)
                        extents = []
                        actual = None if only else path.read_bytes()
                        for rank, _, line in entries:
                            values = dict(item.split("=") for item in line.split())
                            rank, ident = int(rank), int(values["chunk"])
                            chunk = (rank * 20 + 12 - ident) // 7
                            expected = b"".join(records(rank, chunk))
                            begin, size = int(values["offset"]), int(values["bytes"])
                            self.assertEqual(size, len(expected))
                            if actual is not None:
                                self.assertEqual(actual[begin:begin + size], expected)
                            extents.append((begin, begin + size))
                        ordered = sorted(extents)
                        self.assertEqual(ordered[0][0], 0)
                        for left, right in zip(ordered, ordered[1:]):
                            self.assertEqual(left[1], right[0])
                        if actual is not None:
                            self.assertEqual(len(actual), ordered[-1][1])
                        else:
                            self.assertFalse(path.exists())
                        reservations = re.findall(r"^\s+reservations\s+(\d+)\s*$", run.stderr, re.M)
                        self.assertEqual(sorted(map(int, reservations)), [0, 2, 2])
                        if only:
                            samples = re.findall(r"\[SWBWA chunk sample rank (\d+)/(\d+)\] (.*)", run.stderr)
                            self.assertEqual(len(samples), 4)
                            for rank, _, line in samples:
                                values = dict(item.split("=") for item in line.split())
                                rank, ident = int(rank), int(values['chunk'])
                                chunk = (rank * 20 + 12 - ident) // 7
                                blobs = records(rank, chunk)[:100]
                                hashes = [reference_hash(blob) for blob in blobs]
                                xor = 0
                                for value in hashes:
                                    xor ^= value
                                self.assertEqual(int(values['sum'], 16), sum(hashes) & MASK)
                                self.assertEqual(int(values['xor'], 16), xor)

    def test_same_extents_without_file_and_chunk_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            md5_link = build_md5_object(directory)
            outputs = []
            for only in (0, 1):
                exe = directory / f"writer{only}"
                subprocess.run(shlex.split(os.environ.get("MPICC", "mpicc")) + [
                    "-std=gnu11", "-Wall", "-Wextra", "-Werror", "-Wno-unused-function",
                    "-Iinclude", "-DSWBWA_USE_MPI=1", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                    "-DSWBWA_EXEC_MODE=SWBWA_EXEC_SINGLE_CG",
                    "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_UNORDERED",
                    "-DSWBWA_OUTPUT_BUFFER_BYTES=128", f"-DSWBWA_OUTPUT_RMA_ONLY={only}",
                    "tests/test_output_rma.c", "src/host/swbwa_output.c", "-o", str(exe)] + md5_link,
                    cwd=ROOT, check=True)
                path = directory / f"result{only}.sam"
                run = subprocess.run(shlex.split(os.environ.get("MPIEXEC", "mpiexec")) + [
                    "-n", "3", str(exe), str(path)], capture_output=True, text=True, check=True)
                outputs.append(run.stderr)
                if only:
                    self.assertFalse(path.exists())
                else:
                    expected = b"".join(b"".join(records(r, c)) for r in range(2) for c in range(2))
                    actual = path.read_bytes()
                    self.assertEqual(len(actual), len(expected))
                    self.assertEqual(sorted(actual.splitlines()), sorted(expected.splitlines()))
            for label in ("buffered flushes", "buffered flush bytes", "direct oversized writes",
                          "direct oversized bytes", "reservations"):
                pattern = rf"^\s+{label}\s+(\d+)\s*$"
                self.assertEqual(re.findall(pattern, outputs[0], re.M),
                                 re.findall(pattern, outputs[1], re.M), label)
            samples = re.findall(r"\[SWBWA chunk sample rank (\d+)/(\d+)\] (.*)", outputs[1])
            self.assertEqual(len(samples), 4)
            for rank, _, line in samples:
                values = dict(item.split("=") for item in line.split())
                chunk = int(values["chunk"]) % 2
                blobs = records(int(rank), chunk)
                hashes = [reference_hash(blob) for blob in blobs[:100]]
                xor = 0
                for value in hashes:
                    xor ^= value
                self.assertEqual(int(values["sample_reads"]), min(100, len(blobs)))
                self.assertEqual(int(values["sample_bytes"]), sum(map(len, blobs[:100])))
                self.assertEqual(int(values["sam_bytes"]), sum(map(len, blobs)))
                self.assertEqual(int(values["sum"], 16), sum(hashes) & MASK)
                self.assertEqual(int(values["xor"], 16), xor)
            summaries = re.findall(r"\[SWBWA output RMA-only rank (\d+)/(\d+)\] (.*)", outputs[1])
            self.assertEqual(len(summaries), 3)
            for rank, _, line in summaries:
                values = dict(item.split("=") for item in line.split())
                self.assertEqual(values["sam_bytes"], values["reserved_bytes"])
                self.assertEqual(int(values["disk_bytes"]), 0)
                if int(rank) == 2:
                    self.assertEqual(int(values["chunks"]), 0)
                    self.assertEqual(int(values["reservations"]), 0)


if __name__ == "__main__":
    unittest.main()
