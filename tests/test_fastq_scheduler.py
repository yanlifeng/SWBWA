"""Native MPI coverage/reopen tests for both dynamic ticket policies."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SchedulerTest(unittest.TestCase):
    def test_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for exact in (0, 1):
                exe = tmp / f"scheduler-{exact}"
                subprocess.run(shlex.split(os.getenv("MPICC", "mpicc")) + [
                    "-std=gnu11", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-Wno-unused-function", "-Iinclude", "-DSWBWA_USE_MPI=1",
                    "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                    "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_DYNAMIC",
                    f"-DSWBWA_MPI_EXACT_READ_INDEX={exact}",
                    "tests/test_fastq_scheduler.c", "src/host/swbwa_mpi.c",
                    "-pthread", "-o", str(exe)], cwd=ROOT, check=True)
                for mode in ("distributed", "global"):
                    for records, size, tail in ((0, 173, 0), (1, 100000, 0),
                                                (97, 173, 0), (11, 5, 0),
                                                (97, 173, 10)):
                        with self.subTest(exact=exact, mode=mode, records=records,
                                          size=size, tail=tail):
                            env = dict(os.environ, SWBWA_MPI_TICKET_MODE=mode,
                                       SWBWA_MPI_TAIL_PERCENT=str(tail))
                            result = subprocess.run(shlex.split(os.getenv("MPIEXEC", "mpiexec")) + [
                                "-n", "3", str(exe), str(tmp / "input.fq"),
                                str(size), str(records)], env=env, capture_output=True,
                                text=True, timeout=60)
                            self.assertEqual(result.returncode, 0, result.stderr)
                            self.assertIn("coverage PASS", result.stdout)
                            if mode == "global":
                                counts = re.findall(r"tickets beyond the final chunk\s+(\d+)",
                                                    result.stderr)
                                self.assertEqual(list(map(int, counts)), [1] * 6)
                for bad_mode in ("invalid", "mixed"):
                    result = subprocess.run(shlex.split(os.getenv("MPIEXEC", "mpiexec")) + [
                        "-n", "3", str(exe), str(tmp / "input.fq"), "173", "97",
                        bad_mode], capture_output=True, text=True, timeout=60)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("policy rejection PASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
