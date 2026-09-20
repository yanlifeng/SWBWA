"""Local-only checks. Submission tests use fake bsub/bjobs, never SSH."""
import ast
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ReviewToolsTests(unittest.TestCase):
    def test_syntax(self):
        files = [ROOT / "scripts/correctness.sh", ROOT / "scripts/prep_data.sh",
                 ROOT / "scripts/startup_driver3.sh"]
        for directory in ("scripts/ab_profile", "scripts/deploy_ab", "tests"):
            files.extend((ROOT / directory).rglob("*.sh"))
            for path in (ROOT / directory).rglob("*.py"):
                ast.parse(path.read_text(), filename=str(path))
        for path in files:
            subprocess.run(["bash", "-n", str(path)], check=True)

    def test_make_configuration_matrix(self):
        for execution in ("single", "cgs", "cgs_cross"):
            for allocator in ("system", "pool"):
                for mpi in (0, 1):
                    args = ["make", "-s", "print-config", f"EXEC_MODE={execution}",
                            f"CPE_ALLOCATOR={allocator}", f"USE_MPI={mpi}"]
                    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
        for tail in ("0", "10", "100", "", "-1", "101", "010", "1x"):
            result = subprocess.run(
                ["make", "-s", "print-config", "USE_MPI=1", f"MPI_TAIL_PERCENT={tail}"],
                cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode == 0, tail in ("0", "10", "100"))
        for cap in ("0", "64", "", "-1", "064", "x", "2147483648"):
            result = subprocess.run(
                ["make", "-s", "print-config", "USE_MPI=1", f"DISCARD_HASH_BYTES={cap}"],
                cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode == 0, cap in ("0", "64"))

    def test_ab_condition_stops_after_detachment_or_missing_report(self):
        for detached in (True, False):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                binary = root / "runs/old_cross_p0/SWBWA"
                binary.parent.mkdir(parents=True)
                binary.touch(mode=0o755)
                data = root / "data"
                (data / "bwa_test_big_data").mkdir(parents=True)
                (data / "GRCh38.d1.vd1.fa").touch()
                (data / "bwa_test_big_data/ERR1203383_1.fastq").touch()
                stubs = root / "bin"
                stubs.mkdir()
                fake = "#!/bin/bash\necho called >> \"$AB_ROOT/calls\"\necho 'Job <42> has been submitted'\n"
                if not detached:
                    fake += "echo 'Job 42 has been finished.'\n"
                (stubs / "bsub").write_text(fake)
                (stubs / "bsub").chmod(0o755)
                (stubs / "bjobs").write_text("#!/bin/bash\nexit 0\n")
                (stubs / "bjobs").chmod(0o755)
                env = dict(os.environ, AB_ROOT=str(root), DATA_ROOT=str(data),
                           PATH=str(stubs) + os.pathsep + os.environ["PATH"])
                result = subprocess.run(
                    ["bash", str(ROOT / "scripts/ab_profile/ab_cond.sh"), "cross",
                     "ERR1203383", "SE", str(root / "logs")],
                    env=env, capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual((root / "calls").read_text().splitlines(), ["called"])


if __name__ == "__main__":
    unittest.main()
