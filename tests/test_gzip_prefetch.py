#!/usr/bin/env python3
"""Local Source lifecycle regression; no MPI launcher or codec execution."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sanitize", choices=("address,undefined", "thread"))
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="swbwa-prefetch-test-") as directory:
        binary = Path(directory) / "test"
        flags = ["-std=c++11", "-O1", "-g", "-Wall", "-Wextra", "-pthread",
                 "-D_FILE_OFFSET_BITS=64", "-D_POSIX_C_SOURCE=200809L",
                 "-I" + str(ROOT / "include")]
        if args.sanitize:
            flags += ["-fsanitize=" + args.sanitize, "-fno-sanitize-recover=all",
                      "-fno-omit-frame-pointer", "-fno-pie", "-no-pie"]
        subprocess.run([os.getenv("MPICXX", "mpicxx"), *flags,
                        str(ROOT / "tests/test_gzip_prefetch.cpp"),
                        "-Wl,--wrap=pread64", "-Wl,--wrap=MPI_Wtime", "-lz",
                        "-o", str(binary)], check=True, timeout=90)
        subprocess.run([str(binary)], check=True, timeout=30)


if __name__ == "__main__":
    main()
