"""Exercise production B and optional C with the actual allocator under sanitizers."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    flags = shlex.split(os.environ.get("CC", "cc")) + [
        "-std=gnu11", "-O2", "-g", "-Wall", "-Wextra",
        "-Wno-unused-parameter", "-Wno-sign-compare",
        "-fsanitize=address,undefined", "-fno-sanitize-recover=all", "-fno-pie", "-no-pie",
        "-DSWBWA_EXEC_MODE=3", "-DSWBWA_CPE_ALLOC_MODE=2", "-DSWBWA_USE_MPI=0",
        "-include", str(ROOT / "include/swbwa_config.h"),
        "-I", str(ROOT / "src/slave"), "-I", str(ROOT / "include"),
        "-I", str(ROOT / "tests/stubs")]
    with tempfile.TemporaryDirectory(prefix="ldm-policy-check-") as temp:
        executable = str(Path(temp) / "check")
        subprocess.run(flags + [str(ROOT / "tests/test_ldm_policy_allocator.c"),
                               str(ROOT / "src/slave/ldm_alloc.c"),
                               str(ROOT / "src/slave/malloc_wrap.c"),
                               "-lm", "-o", executable], check=True)
        # No mode macro or policy arguments: verify the real production defaults.
        subprocess.run([executable], check=True)
        subprocess.run([executable, "0,4096,1024,0,0,0,0,0,0,0,0",
                        "0,4096,4096,0,4096,0,0,0,0,0,0"], check=True)


if __name__ == "__main__":
    main()
