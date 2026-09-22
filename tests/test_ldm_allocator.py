"""ASan/UBSan allocator and actual ksw_global2 differential tests, no Sunway SDK."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = (ROOT / "src/slave/ksw.c").read_text()
    start = source.index("static inline uint32_t *push_cigar(")
    end = source.index("\nint ksw_global(", start)
    prefix = """
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include "malloc_wrap.h"
typedef struct { int32_t h, e; } eh_t;
typedef struct { int32_t lane[16]; } __m128i;
typedef struct _kswq_t kswq_t;
#define SWBWA_KSW_U8_LANES 16
#define SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES 16384
#define MINUS_INF (-0x40000000)
#define LIKELY(x) (x)
void *swbwa_test_unknown_malloc(int *bytes) { return malloc((*bytes)++); }
void *swbwa_test_unknown_calloc(int *bytes) { return calloc(1, (*bytes)++); }
char *swbwa_test_unknown_strdup(const char *s) { return strdup(s); }
"""
    canonical = None
    with tempfile.TemporaryDirectory(prefix="swbwa-ldm-allocator-") as tmp:
        tmp = Path(tmp)
        kernel = tmp / "global.c"
        query_struct = source[source.index("struct _kswq_t {"):source.index("\n#ifndef SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES")]
        query_code = source[source.index("static kswq_t *ksw_qinit_impl("):source.index("\n#if defined __ARM_NEON", source.index("static void ksw_qdestroy"))]
        query_test = """
void *swbwa_test_cached_query(void) {
    uint8_t query[8] = {0,1,2,3,0,1,2,3};
    int8_t mat[25];
    for (int i = 0; i < 25; ++i) mat[i] = i / 5 == i % 5 ? 1 : -4;
    return ksw_qinit(1, 8, query, 5, mat);
}
void swbwa_test_cached_query_free(void *p) {
    kswq_t *q = p;
    assert(q->qlen == 8 && q->in_ldm == 0);
    assert(((int *)q->qp)[0] == 5);
    ksw_qdestroy(q);
}
"""
        kernel.write_text(prefix + query_struct + query_code + query_test + source[start:end])
        for mode, cap, manual, budget in (
            (1, 4096, 1, 40960), (2, 4096, 1, 40960),
            (2, 1024, 1, 40960), (2, 256, 1, 40960),
            (2, 4480, 0, 40960), (2, 32768, 0, 40960),
            (2, 98304, 0, 114688),
            (3, 256, 0, 40960), (3, 4096, 1, 40960),
            (3, 32768, 0, 40960), (3, 98304, 0, 114688),
            (4, 256, 0, 40960), (4, 4096, 1, 40960),
            (4, 32768, 0, 40960), (4, 98304, 0, 114688),
            (2, 204800, 0, 229376), (3, 204800, 0, 229376),
            (4, 204800, 0, 229376),
            (4, 188416, 0, 212992),
        ):
            exe = tmp / f"test-{mode}-{cap}"
            # Keep regression coverage for historical allocator policies without
            # exposing their old combinations as production build switches.
            config = tmp / "allocator-test-config.h"
            config.write_text(
                '#define SWBWA_CPE_LDM_MODE 1\n'
                '#include "swbwa_config.h"\n'
                '#undef SWBWA_CPE_LDM_ALLOC\n'
                '#undef SWBWA_CPE_MANUAL_LDM\n'
                f'#define SWBWA_CPE_LDM_ALLOC {mode}\n'
                f'#define SWBWA_CPE_MANUAL_LDM {manual}\n')
            cmd = shlex.split(os.environ.get("CC", "cc")) + [
                "-std=gnu11", "-O2", "-g", "-Wall", "-Wextra",
                "-Wno-unused-parameter", "-Wno-sign-compare",
                "-Werror=implicit-function-declaration", "-fsanitize=address,undefined",
                "-fno-sanitize-recover=all", "-fno-pie", "-no-pie",
                "-DSWBWA_EXEC_MODE=3", "-DSWBWA_CPE_ALLOC_MODE=2",
                f"-DSWBWA_CPE_LDM_BYTES={cap}",
                f"-DSWBWA_LDM_SCRATCH_BUDGET_BYTES={budget}",
                "-include", str(config),
                "-I", str(ROOT / "src/slave"), "-I", str(ROOT / "include"),
                "-I", str(ROOT / "tests/stubs"),
                str(ROOT / "tests/test_ldm_allocator.c"), str(kernel),
                str(ROOT / "src/slave/ldm_alloc.c"),
                str(ROOT / "src/slave/malloc_wrap.c"), "-lm", "-o", str(exe)]
            subprocess.run(cmd, check=True)
            result = subprocess.check_output([str(exe)])
            if canonical is None:
                canonical = result
            assert canonical == result, (mode, cap, "DP scores/CIGAR differ")
            print(f"mode={mode} capacity={cap} manual={manual} budget={budget}: "
                  f"{result.decode().splitlines()[-1]}", flush=True)


if __name__ == "__main__":
    main()
