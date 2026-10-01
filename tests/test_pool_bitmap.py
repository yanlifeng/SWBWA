"""Reproduce the final candidate's bitmap and epoch safety checks locally."""
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
        "-Werror=implicit-function-declaration", "-fsanitize=address,undefined",
        "-fno-sanitize-recover=all", "-fno-pie", "-no-pie",
        "-DSWBWA_EXEC_MODE=3", "-DSWBWA_CPE_ALLOC_MODE=2",
        "-DSWBWA_CPE_LDM_MODE=3", "-include", str(ROOT / "include/swbwa_config.h"),
        "-I", str(ROOT / "src/slave"), "-I", str(ROOT / "include"),
        "-I", str(ROOT / "tests/stubs"),
    ]
    with tempfile.TemporaryDirectory(prefix="swbwa-bitmap-check-") as name:
        tmp = Path(name)
        differential = tmp / "differential"
        subprocess.run(flags + [str(ROOT / "tests/test_pool_bitmap.c"),
                       str(ROOT / "src/slave/ldm_alloc.c"), "-lm", "-o",
                       str(differential)], check=True)
        subprocess.run([str(differential)], check=True)
        epoch = tmp / "epoch"
        source = (ROOT / "tests/test_ldm_heap_cache.c").read_text()
        marker = "    assert(!backing); backing = malloc(bytes);"
        assert source.count(marker) == 1
        source = source.replace(marker, '    if (getenv("TEST_REFUSE_LDM")) return NULL;\n' + marker)
        harness = tmp / "epoch.c"
        harness.write_text(source)
        subprocess.run(flags + [str(harness),
                       str(ROOT / "src/slave/malloc_wrap.c"),
                       str(ROOT / "src/slave/ldm_alloc.c"), "-lm", "-o",
                       str(epoch)], check=True)
        for mode in (0, 1, 2):
            print(f"heap epoch test: bitmap mode {mode}", flush=True)
            subprocess.run([str(epoch)], check=True,
                           env=dict(os.environ, TEST_POOL_BITMAP=str(mode)))
        print("bitmap mode 2 with complete SDK LDM refusal", flush=True)
        subprocess.run([str(epoch)], check=True,
                       env=dict(os.environ, TEST_POOL_BITMAP="2", TEST_REFUSE_LDM="1"))
        tree = tmp / "chain_tree"
        subprocess.run(flags + [str(ROOT / "tests/test_ldm_chain_tree.c"),
                       str(ROOT / "src/slave/malloc_wrap.c"),
                       str(ROOT / "src/slave/ldm_alloc.c"), "-lm", "-o", str(tree)],
                       check=True)
        subprocess.run([str(tree)], check=True)


if __name__ == "__main__":
    main()
