"""Check Make defaults without invoking the Sunway compiler or linker."""
import itertools
import os
from pathlib import Path
import shlex
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OptimizationDefaultsTest(unittest.TestCase):
    def config(self, **options):
        env = {key: value for key, value in os.environ.items()
               if key not in {"MAKEFLAGS", "MFLAGS", "MAKEOVERRIDES",
                              "CPE_KERNEL_OPT", "CPE_DISCARD_DIGEST",
                              "CPE_PROFILE", "HOST_MALLOC_STATS", "HOST_MPE_THREADS",
                              "CPE_LDM_MODE"}}
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
                self.assertEqual(values["CPE_LDM_MODE"], "3" if core else "2")
                self.assertEqual(values["CPE_DISCARD_DIGEST"], str(int(digest)))
                self.assertEqual(values["HOST_MPE_THREADS"], "1" if execution == "single" else "6")

    def test_host_workers_override(self):
        for execution in ("single", "cgs", "cgs_cross"):
            for workers in (1, 6):
                with self.subTest(execution=execution, workers=workers):
                    result = self.config(EXEC_MODE=execution, HOST_MPE_THREADS=workers)
                    if execution == "single" and workers == 6:
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("HOST_MPE_THREADS=6 requires", result.stderr)
                    else:
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertIn(f"HOST_MPE_THREADS={workers}\n", result.stdout)

    def test_header_host_workers_defaults_and_overrides(self):
        cc = shlex.split(os.environ.get("CC", "cc"))
        for execution, mpi, workers in itertools.product((1, 2, 3), (0, 1), (None, 1, 6)):
            with self.subTest(execution=execution, mpi=mpi, workers=workers):
                flags = [f"-DSWBWA_EXEC_MODE={execution}", f"-DSWBWA_USE_MPI={mpi}"]
                if workers is not None:
                    flags += [f"-DSWBWA_HOST_MPE_THREADS={workers}"]
                expected = workers if workers is not None else (1 if execution == 1 else 6)
                source = ('#include "swbwa_config.h"\n'
                          f'_Static_assert(SWBWA_HOST_MPE_THREADS == {expected}, "MPE count");\n')
                result = subprocess.run(cc + ["-std=c11", "-Iinclude", *flags,
                                              "-fsyntax-only", "-x", "c", "-"],
                                        cwd=ROOT, input=source, text=True, capture_output=True)
                if execution == 1 and workers == 6:
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("SWBWA_HOST_MPE_THREADS=6 requires", result.stderr)
                else:
                    self.assertEqual(result.returncode, 0, result.stderr)

    def test_explicit_disable(self):
        result = self.config(EXEC_MODE="cgs_cross", CPE_ALLOCATOR="pool",
                             USE_MPI=0, OUTPUT_MODE="discard", DISCARD_HASH_BYTES=0,
                             CPE_KERNEL_OPT=0, CPE_DISCARD_DIGEST=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("CPE_KERNEL_OPT=0\n", result.stdout)
        self.assertIn("CPE_DISCARD_DIGEST=0\n", result.stdout)

    def test_retired_cpe_experiments_are_rejected(self):
        cc = shlex.split(os.environ.get("CC", "cc"))
        for option in ("SWBWA_CPE_POOL_INDEX", "SWBWA_CPE_POOL_INDEX_CAPACITY",
                       "SWBWA_CPE_TASK_LEASE", "SWBWA_CPE_TASK_QUEUES",
                       "SWBWA_CROSS_FINISH_BARRIER", "SWBWA_CROSS_POLL_BACKOFF",
                       "SWBWA_ENABLE_PACKED_INT8"):
            with self.subTest(option=option):
                result = subprocess.run(
                    cc + ["-std=c11", "-Iinclude", f"-D{option}=1",
                          "-fsyntax-only", "-x", "c", "-"],
                    cwd=ROOT, input='#include "swbwa_config.h"\n',
                    text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("retired CPE experiment option", result.stderr)

    def test_build_layout(self):
        for execution, mpi in itertools.product(("single", "cgs", "cgs_cross"), (0, 1)):
            with self.subTest(execution=execution, mpi=mpi):
                result = subprocess.run(
                    ["make", "-Bn", "--no-print-directory", "SWBWA",
                     f"EXEC_MODE={execution}", f"USE_MPI={mpi}",
                     "CPE_ALLOCATOR=pool", "OUTPUT_MODE=split"],
                    cwd=ROOT, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("libswbwa.a", result.stdout)
                self.assertIn("-lswbwa", result.stdout)
                self.assertNotIn("-lbwa ", result.stdout)
                for removed in ("src/slave/main.o", "src/slave/example.o",
                                "src/slave/kthread.o", "src/host/example.o"):
                    self.assertNotIn(removed, result.stdout)

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
