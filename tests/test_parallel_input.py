"""Byte-exact positioned input, partial transfers, EINTR and worker bounds."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ParallelInputTest(unittest.TestCase):
    def test_positioned_reads(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            for execution, workers in ((1, None), (2, None), (3, None), (3, 1), (3, 6)):
                for wrapped in (False, True):
                    with self.subTest(execution=execution, workers=workers, wrapped=wrapped):
                        exe = directory / f"read-{execution}-{workers}-{int(wrapped)}"
                        args = shlex.split(os.environ.get("CC", "cc")) + [
                            "-std=gnu11", "-O2", "-g", "-D_GNU_SOURCE", "-D_FILE_OFFSET_BITS=64",
                            "-Wall", "-Wextra", "-Werror", "-fsanitize=address,undefined",
                            "-fno-sanitize-recover=all",
                            "-fno-omit-frame-pointer", "-fno-pie", "-no-pie",
                            f"-DSWBWA_EXEC_MODE={execution}", "-Iinclude",
                            "tests/test_parallel_input.c", "src/host/swbwa_input.c",
                            "src/host/swbwa_host_workers.c", "-pthread", "-o", str(exe)]
                        if workers is not None:
                            args += [f"-DSWBWA_HOST_MPE_THREADS={workers}"]
                        if wrapped:
                            args += ["-DTEST_WRAP_READ", "-Wl,--wrap=pread64"]
                            # _FILE_OFFSET_BITS=64 maps calls and wrapper names alike.
                            args += ["-D__wrap_pread=__wrap_pread64", "-D__real_pread=__real_pread64"]
                        subprocess.run(args, cwd=ROOT, check=True)
                        subprocess.run([str(exe), str(directory / "input.bin")], check=True)


if __name__ == "__main__":
    unittest.main()
