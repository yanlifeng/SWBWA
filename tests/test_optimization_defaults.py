"""Check Make defaults without invoking the Sunway compiler or linker."""
import itertools
import os
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OptimizationDefaultsTest(unittest.TestCase):
    def config(self, **options):
        env = {key: value for key, value in os.environ.items()
               if key not in {"MAKEFLAGS", "MFLAGS", "MAKEOVERRIDES",
                              "CPE_KERNEL_OPT", "CPE_DISCARD_DIGEST",
                              "CPE_PROFILE", "HOST_MALLOC_STATS"}}
        return subprocess.run(
            ["make", "-s", "--no-print-directory", "print-config",
             *[f"{key}={value}" for key, value in options.items()]],
            cwd=ROOT, env=env, text=True, capture_output=True)

    def test_default_matrix(self):
        for execution, allocator, mpi, output, prefix in itertools.product(
                ("single", "cgs", "cgs_cross"), ("system", "pool"),
                (0, 1), ("split", "single_unordered", "discard"), (0, 64)):
            with self.subTest(execution=execution, allocator=allocator,
                              mpi=mpi, output=output, prefix=prefix):
                result = self.config(EXEC_MODE=execution, CPE_ALLOCATOR=allocator,
                                     USE_MPI=mpi, OUTPUT_MODE=output,
                                     DISCARD_HASH_BYTES=prefix)
                if not mpi and output == "single_unordered":
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("requires USE_MPI=1", result.stderr)
                    continue
                self.assertEqual(result.returncode, 0, result.stderr)
                values = dict(line.split("=", 1) for line in result.stdout.splitlines())
                core = execution == "cgs_cross" and allocator == "pool" and mpi == 0
                digest = core and output == "discard" and prefix == 0
                self.assertEqual(values["CPE_KERNEL_OPT"], str(int(core)))
                self.assertEqual(values["CPE_DISCARD_DIGEST"], str(int(digest)))

    def test_explicit_disable(self):
        result = self.config(EXEC_MODE="cgs_cross", CPE_ALLOCATOR="pool",
                             USE_MPI=0, OUTPUT_MODE="discard", DISCARD_HASH_BYTES=0,
                             CPE_KERNEL_OPT=0, CPE_DISCARD_DIGEST=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("CPE_KERNEL_OPT=0\n", result.stdout)
        self.assertIn("CPE_DISCARD_DIGEST=0\n", result.stdout)

    def test_explicit_digest_rejects_unsupported_modes(self):
        target = dict(EXEC_MODE="cgs_cross", CPE_ALLOCATOR="pool", USE_MPI=0,
                      OUTPUT_MODE="discard", DISCARD_HASH_BYTES=0,
                      CPE_DISCARD_DIGEST=1)
        for override in ({"EXEC_MODE": "single"}, {"EXEC_MODE": "cgs"},
                         {"CPE_ALLOCATOR": "system"}, {"USE_MPI": 1},
                         {"DISCARD_HASH_BYTES": 64}):
            with self.subTest(override=override):
                result = self.config(**dict(target, **override))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("CPE_DISCARD_DIGEST=1 requires", result.stderr)

    def test_invalid_boolean_values(self):
        for key, value in itertools.product(
                ("CPE_KERNEL_OPT", "CPE_DISCARD_DIGEST"),
                ("", "2", "-1", "0 1", "1 junk", "%")):
            with self.subTest(key=key, value=value):
                result = self.config(EXEC_MODE="cgs_cross", CPE_ALLOCATOR="pool",
                                     USE_MPI=0, OUTPUT_MODE="discard",
                                     DISCARD_HASH_BYTES=0, **{key: value})
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"{key} must be 0 or 1", result.stderr)


if __name__ == "__main__":
    unittest.main()
