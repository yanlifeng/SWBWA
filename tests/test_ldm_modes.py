"""Check the three production modes and the tracked SDK allocation gate."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CC = shlex.split(os.environ.get("CC", "cc"))
BASE = ["-DSWBWA_EXEC_MODE=3", "-DSWBWA_CPE_ALLOC_MODE=2", "-DSWBWA_USE_MPI=0"]


def main():
    for mode in range(3):
        result = subprocess.check_output([
            "make", "-s", "print-config", "EXEC_MODE=cgs_cross", "CPE_ALLOCATOR=pool",
            "USE_MPI=0", f"CPE_LDM_MODE={mode}"], cwd=ROOT, text=True)
        assert f"CPE_LDM_MODE={mode}" in result
    for bad in ("CPE_LDM_MODE=3", "CPE_LDM_MODE=", "CPE_LDM_MODE=0 1", "CPE_LDM_MODE=%",
                "CPE_LDM_ALLOC=off", "CPE_MANUAL_LDM=0", "CPE_LDM_BYTES=32768",
                "LDM_SCRATCH_BUDGET=40960"):
        assert subprocess.run(["make", "-s", "print-config", bad], cwd=ROOT,
                              capture_output=True).returncode != 0, bad

    with tempfile.TemporaryDirectory(prefix="swbwa-ldm-modes-") as tmp:
        tmp = Path(tmp)
        test = tmp / "gate.c"
        test.write_text(r'''
#include <assert.h>
#include <stdlib.h>
int swbwa_test_cpe_id;
static int allocations;
void *ldm_malloc(unsigned long n) { ++allocations; return malloc(n); }
void ldm_free(void *p, unsigned long n) { (void)n; free(p); }
void swbwa_cpe_fail(int c, long a, long b, long d) { abort(); }
int main(void) {
    assert(SWBWA_CPE_LDM_ALLOC == (SWBWA_CPE_LDM_MODE == 1 ? 4 : 0));
    assert(SWBWA_CPE_MANUAL_LDM == (SWBWA_CPE_LDM_MODE == 2));
    assert(SWBWA_CPE_LDM_BYTES == 32768);
    swbwa_ldm_begin_batch();
    for (int site = 0; site <= 9; ++site) {
        void *p = swbwa_ldm_alloc(64, site);
        int expected = SWBWA_CPE_LDM_MODE == 2 ||
                      (SWBWA_CPE_LDM_MODE == 1 && site == SWBWA_LDM_AUTO_ARENA_SITE);
        assert((p != NULL) == expected);
        swbwa_ldm_release(p, 64);
    }
    assert(allocations == (SWBWA_CPE_LDM_MODE == 0 ? 0 :
                          SWBWA_CPE_LDM_MODE == 1 ? 1 : 10));
    assert(swbwa_ldm_outstanding() == 0);
    return 0;
}
''')
        common = CC + ["-std=gnu11", "-O2", "-g", "-fsanitize=address,undefined",
                       "-fno-sanitize-recover=all", "-fno-pie", "-no-pie"] + BASE + [
            "-include", str(ROOT / "include/swbwa_config.h"),
            "-I", str(ROOT / "include"), "-I", str(ROOT / "tests/stubs")]
        for mode in range(3):
            exe = tmp / f"mode{mode}"
            subprocess.run(common + [f"-DSWBWA_CPE_LDM_MODE={mode}", str(test),
                           str(ROOT / "src/slave/malloc_wrap.c"), "-lm", "-o", str(exe)], check=True)
            subprocess.run([str(exe)], check=True)
            print(f"CPE_LDM_MODE={mode}: PASS configuration + SDK gate", flush=True)
        for flags in (["-DSWBWA_CPE_LDM_MODE=3"], ["-DSWBWA_CPE_LDM_ALLOC=0"],
                      ["-DSWBWA_CPE_MANUAL_LDM=0"],
                      ["-DSWBWA_CPE_LDM_MODE=1", "-USWBWA_USE_MPI", "-DSWBWA_USE_MPI=1"]):
            assert subprocess.run(CC + BASE + flags + ["-fsyntax-only", "-x", "c",
                   "-include", str(ROOT / "include/swbwa_config.h"), "-"],
                   input="", text=True, capture_output=True).returncode != 0, flags
    print("PASS retired-option and incompatible-mode checks")


if __name__ == "__main__":
    main()
