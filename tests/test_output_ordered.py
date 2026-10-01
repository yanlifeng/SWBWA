"""Real MPI prefix ordering, empty chunks, repeated open and no-file mode."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest
from test_sam_stream_md5 import build_md5_object

ROOT = Path(__file__).resolve().parents[1]


class OrderedOutputTest(unittest.TestCase):
    def test_prefix_and_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            md5_link = build_md5_object(tmp)
            for exact, only, workers in ((0, 0, 1), (1, 0, 6), (0, 1, 6)):
                exe = tmp / f"ordered-{exact}-{only}-{workers}"
                subprocess.run(shlex.split(os.getenv("MPICC", "mpicc")) + [
                    "-std=gnu11", "-O2", "-Wall", "-Wextra", "-Werror", "-Wno-unused-function",
                    "-Iinclude", "-DSWBWA_USE_MPI=1", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                    "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_DYNAMIC",
                    "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_ORDERED",
                    "-DSWBWA_OUTPUT_BUFFER_BYTES=16", f"-DSWBWA_OUTPUT_RMA_ONLY={only}",
                    f"-DSWBWA_MPI_EXACT_READ_INDEX={exact}", f"-DSWBWA_HOST_MPE_THREADS={workers}",
                    "tests/test_output_ordered.c", "src/host/swbwa_mpi.c",
                    "src/host/swbwa_output.c", "src/host/swbwa_host_workers.c", "-pthread",
                    "-o", str(exe)] + md5_link, cwd=ROOT, check=True)
                for records, size, tail in ((0, 173, 0), (1, 10000, 0),
                                            (97, 173, 0), (11, 5, 0), (97, 173, 10)):
                    with self.subTest(exact=exact, only=only, records=records, size=size, tail=tail):
                        path = tmp / f"out-{exact}-{only}-{records}-{size}-{tail}.sam"
                        env = dict(os.environ, SWBWA_MPI_TAIL_PERCENT=str(tail))
                        env.pop("SWBWA_MPI_TICKET_MODE", None)
                        run = subprocess.run(shlex.split(os.getenv("MPIEXEC", "mpiexec")) + [
                            "-n", "3", str(exe), str(tmp / "input.fq"), str(path), str(size),
                            str(records)], capture_output=True, text=True, env=env, timeout=60)
                        self.assertEqual(run.returncode, 0, run.stderr)
                        self.assertIn("ordered output PASS", run.stdout)
                        expected = b"".join(f"read{i:03d}\t0\tchr1\t{i+1}\nread{i:03d}\t2048\tchr2\t{i+2}\n".encode()
                                            for i in range(records))
                        if only:
                            self.assertFalse(path.exists())
                        else:
                            self.assertEqual(path.read_bytes(), expected)
                        extents = re.findall(r"\[SWBWA output extent rank \d+/\d+\] (.*)", run.stderr)
                        entries = [dict((k, int(v)) for k, v in (word.split("=") for word in line.split()))
                                   for line in extents]
                        unique = {entry['chunk']: entry for entry in entries}
                        # Two opens truncate and reproduce the same extent map.
                        self.assertEqual(len(entries), 2 * len(unique))
                        offset = 0
                        for entry in sorted(unique.values(), key=lambda e: e['chunk']):
                            self.assertEqual(entry['offset'], offset)
                            offset += entry['bytes']
                        self.assertEqual(offset, len(expected))


if __name__ == "__main__":
    unittest.main()
