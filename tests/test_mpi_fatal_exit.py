#!/usr/bin/env python3
"""Exercise SWBWA's MPI startup and fatal exits locally, without a CPE runtime."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    reports = []
    for sunway_guard in (False, True):
        binary = work / ("fatal-sunway-guard" if sunway_guard else "fatal-portable")
        flags = ["-O2", "-g", "-Wall", "-Wextra", "-D_GNU_SOURCE", "-D_FILE_OFFSET_BITS=64",
                 "-ffunction-sections", "-fdata-sections", "-DSWBWA_USE_MPI=1",
                 "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_STATIC",
                 "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SPLIT", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0"]
        if sunway_guard:
            flags.append("-D__sw_host__")
        subprocess.run([os.environ.get("MPICC", "mpicc"), *flags,
                        "-I" + str(ROOT / "include"), str(ROOT / "src/host/swbwa_mpi.c"),
                        str(ROOT / "tests/test_mpi_fatal_exit.c"),
                        "-Wl,--gc-sections", "-pthread", "-o", str(binary)], check=True)
        for mode in ("normal", "world", "dup", "all", "thread"):
            env = dict(os.environ, SWCH_HANG_WHEN_ERROR="1")
            begin = time.monotonic()
            run = subprocess.run(["mpirun", "-launcher", "fork", "-n", "3", str(binary), mode],
                                 env=env, capture_output=True, text=True, timeout=30)
            elapsed = time.monotonic() - begin
            record = dict(sunway_guard=sunway_guard, mode=mode, returncode=run.returncode,
                          elapsed_s=elapsed, stdout=run.stdout, stderr=run.stderr)
            reports.append(record)
            (work / "results.json").write_text(json.dumps(reports, indent=2) + "\n")
            assert (run.returncode == 0) == (mode == "normal"), record
            assert "unexpected return" not in run.stderr, record
            if mode == "normal":
                assert "normal shutdown" in run.stdout, record
            else:
                assert "entering MPI_Abort" in run.stderr, record
            print("PASS", sunway_guard, mode, f"{elapsed:.3f}s", flush=True)


if __name__ == "__main__":
    main()
