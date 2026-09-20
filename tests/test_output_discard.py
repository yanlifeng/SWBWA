"""Compile the real non-MPI discard writer and check its full-byte fingerprint."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MASK = (1 << 64) - 1


def reference_hash(data):
    multiplier = 0xC6A4A7935BD1E995
    value = 0x9E3779B97F4A7C15 ^ ((len(data) * multiplier) & MASK)
    end = len(data) // 8 * 8
    for offset in range(0, end, 8):
        word = int.from_bytes(data[offset:offset + 8], "little")
        word = word * multiplier & MASK
        word ^= word >> 47
        word = word * multiplier & MASK
        value = (value ^ word) * multiplier & MASK
    if end != len(data):
        value ^= int.from_bytes(data[end:], "little")
        value = value * multiplier & MASK
    value ^= value >> 47
    value = value * multiplier & MASK
    return value ^ (value >> 47)


class OutputDiscardTest(unittest.TestCase):
    def test_non_mpi_full_hash_and_no_file(self):
        records = [b"read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"]
        records.extend(bytes(i & 255 for i in range(1, size + 1))
                       for size in range(1, 1025))
        hashes = list(map(reference_hash, records))
        expected_xor = 0
        for value in hashes:
            expected_xor ^= value
        with tempfile.TemporaryDirectory() as directory:
            executable = str(Path(directory) / "discard_test")
            subprocess.run(shlex.split(os.environ.get("CC", "cc")) + [
                "-std=gnu11", "-Wall", "-Wextra", "-Wno-unused-function",
                "-Iinclude", "-DSWBWA_USE_MPI=0",
                "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_DISCARD",
                "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                "tests/test_output_discard.c", "src/host/swbwa_output.c",
                "-o", executable], cwd=ROOT, check=True)
            result = subprocess.run([executable, str(Path(directory) / "out.sam")],
                                    env=dict(os.environ, SWBWA_DISCARD_HASH="1"),
                                    capture_output=True, text=True, check=True)
            report = re.search(r"calls=(\d+) bytes=(\d+) sum=0x([0-9a-f]+) "
                               r"xor=0x([0-9a-f]+) enabled=1 hash_prefix_bytes=0",
                               result.stderr)
            self.assertIsNotNone(report, result.stderr)
            self.assertEqual(int(report[1]), len(records))
            self.assertEqual(int(report[2]), sum(map(len, records)))
            self.assertEqual(int(report[3], 16), sum(hashes) & MASK)
            self.assertEqual(int(report[4], 16), expected_xor)


if __name__ == "__main__":
    unittest.main()
