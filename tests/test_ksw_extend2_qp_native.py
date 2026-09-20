"""Differential host tests of the actual CPE ksw_extend2, QP toggle off/on.

The harness includes src/slave/ksw.c unchanged. Unused SIMD kernels are linked
out; the scalar extension function, scratch storage and QP construction run
natively. This is a correctness test, not a Sunway performance measurement.
LDM placement is left at its production default (zero); the profile assertions
inspect persistent heap scratch and deliberately do not test the LDM candidate.
Run: python3 tests/test_ksw_extend2_qp_native.py
"""

import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
HARNESS = r'''
#include <stdio.h>
#include "src/slave/ksw.c"

_Static_assert(SWBWA_KSW_EXTEND2_LDM_MAX_BYTES == 0,
               "QP-only harness requires the production LDM default to be zero");

static uint32_t state = 7963242;
static unsigned long cases, qp_bytes;

static uint32_t random_word(void)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static void check_case(int qlen, const uint8_t *query, int tlen,
        const uint8_t *target, int m, const int8_t *mat,
        int od, int ed, int oi, int ei, int w, int bonus, int zdrop,
        int h0, unsigned null_mask)
{
    int out[5] = {-777, -777, -777, -777, -777};
    int *ptr[5];
    int score;
    uint8_t saved_query[600], saved_target[600];
    int8_t saved_mat[400];
    size_t live_qp = (size_t)qlen * m;
    size_t live_eh = ((size_t)qlen + 1) * sizeof(eh_t);

    assert(qlen > 0 && qlen <= 600 && tlen >= 0 && tlen <= 600);
    assert(m > 0 && m <= 20 && ed > 0 && ei > 0);
    memcpy(saved_query, query, qlen);
    memcpy(saved_target, target, tlen);
    memcpy(saved_mat, mat, (size_t)m * m);
    if (swbwa_extend2_scratch.qp)
        memset(swbwa_extend2_scratch.qp, 0xa5, swbwa_extend2_scratch.qp_capacity);
    if (swbwa_extend2_scratch.eh)
        memset(swbwa_extend2_scratch.eh, 0x5a, swbwa_extend2_scratch.eh_capacity);
    for (int i = 0; i < 5; ++i) ptr[i] = null_mask & (1u << i) ? NULL : &out[i];
    score = ksw_extend2(qlen, query, tlen, target, m, mat,
                       od, ed, oi, ei, w, bonus, zdrop, h0,
                       ptr[0], ptr[1], ptr[2], ptr[3], ptr[4]);
    ++cases;
    for (int row = 0; row < m; ++row) {
        for (int j = 0; j < qlen; ++j) {
            assert(swbwa_extend2_scratch.qp[(size_t)row * qlen + j]
                   == mat[row * m + query[j]]);
            ++qp_bytes;
        }
    }
    for (size_t i = live_qp; i < swbwa_extend2_scratch.qp_capacity; ++i)
        assert((uint8_t)swbwa_extend2_scratch.qp[i] == 0xa5);
    for (size_t i = live_eh; i < swbwa_extend2_scratch.eh_capacity; ++i)
        assert(((uint8_t *)swbwa_extend2_scratch.eh)[i] == 0x5a);
    for (int i = 0; i < 5; ++i)
        if (null_mask & (1u << i)) assert(out[i] == -777);
    assert(memcmp(query, saved_query, qlen) == 0);
    assert(memcmp(target, saved_target, tlen) == 0);
    assert(memcmp(mat, saved_mat, (size_t)m * m) == 0);
    printf("%lu %d %d %d %d %d %d\n", cases, score,
           out[0], out[1], out[2], out[3], out[4]);
}

static void decode(int value, int len, uint8_t *seq)
{
    for (int i = 0; i < len; ++i) { seq[i] = value % 5; value /= 5; }
}

static void test_short_exhaustive(void)
{
    uint8_t query[3], target[3];
    int8_t mat[25];
    const int powers[] = {1, 5, 25, 125};
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            mat[i * 5 + j] = i == 4 || j == 4 ? -1 : i == j ? 1 : -4;
    for (int qlen = 1; qlen <= 3; ++qlen) {
        for (int q = 0; q < powers[qlen]; ++q) {
            decode(q, qlen, query);
            for (int tlen = 0; tlen <= 3; ++tlen) {
                for (int t = 0; t < powers[tlen]; ++t) {
                    decode(t, tlen, target);
                    check_case(qlen, query, tlen, target, 5, mat,
                               6, 1, 6, 1, 100, 5, 100, 19, 0);
                }
            }
        }
    }
}

static void test_lengths_matrices_and_fallback(void)
{
    const int lengths[] = {1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 60, 62,
                           63, 64, 75, 127, 128, 129, 149, 150, 151, 255, 256, 511};
    const int alphabets[] = {1, 2, 4, 5, 6, 20};
    uint8_t query[600], target[600];
    int8_t mat[400];

    for (int n = 0; n < 4096; ++n) {
        int m = n % 3 ? 5 : alphabets[(n / 3) % 6];
        int qlen = lengths[random_word() % (sizeof(lengths) / sizeof(*lengths))];
        int tlen = n % 17 == 0 ? 0 : random_word() % 301;
        for (int i = 0; i < qlen; ++i) query[i] = random_word() % m;
        for (int i = 0; i < tlen; ++i) target[i] = random_word() % m;
        for (int i = 0; i < m * m; ++i) {
            switch (n % 5) {
            case 0: mat[i] = 0; break; // pervasive ties
            case 1: mat[i] = -128; break;
            case 2: mat[i] = 127; break;
            default: mat[i] = (int)(random_word() % 256) - 128; break;
            }
        }
        int od = random_word() % 20, oi = random_word() % 20;
        int ed = 1 + random_word() % 5, ei = 1 + random_word() % 5;
        int w = random_word() % 201, bonus = (int)(random_word() % 61) - 20;
        int zdrop = n % 4 == 0 ? 0 : 1 + random_word() % 100;
        int h0 = 1 + random_word() % 500;
        check_case(qlen, query, tlen, target, m, mat,
                   od, ed, oi, ei, w, bonus, zdrop, h0, 0);
        if (n % 64 == 0) {
            /* Reuse the same pointers after changing contents, and retry
             * with a different band. Nothing may cache by pointer alone.
             */
            query[0] = (query[0] + 1) % m;
            mat[0] = mat[0] == 127 ? -128 : mat[0] + 1;
            for (unsigned mask = 0; mask < 32; ++mask)
                check_case(qlen, query, tlen, target, m, mat,
                           od, ed, oi, ei, w * 2, bonus, zdrop, h0, mask);
        }
    }
}

int main(void)
{
    test_short_exhaustive();
    test_lengths_matrices_and_fallback();
    fprintf(stderr, "QP_UNROLL=%d: %lu actual extend2 calls, %lu QP bytes checked\n",
            SWBWA_KSW_EXTEND2_QP_UNROLL, cases, qp_bytes);
    free(swbwa_extend2_scratch.qp);
    free(swbwa_extend2_scratch.eh);
    return 0;
}
'''


def main():
    compiler = shlex.split(os.environ.get("CC", "gcc"))
    with tempfile.TemporaryDirectory(prefix="swbwa-extend2-qp-") as temporary:
        path = Path(temporary)
        source = path / "test.c"
        source.write_text(HARNESS, encoding="ascii")
        canonical = None
        for optimization in ("-O2", "-O3"):
            for toggle in (0, 1):
                executable = path / f"qp-{optimization[1:]}-{toggle}"
                command = compiler + [
                    "-std=gnu11", optimization, "-g", "-Wall", "-Wextra",
                    "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections",
                    "-fsanitize=address,undefined", "-fno-sanitize-recover=all",
                    "-fno-omit-frame-pointer", "-fno-pie", "-no-pie",
                    "-DSWBWA_ENABLE_CPE_PROFILE=0",
                    "-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=0",
                    "-DSWBWA_ENABLE_MATESW_DUAL_FORWARD=0",
                    f"-DSWBWA_KSW_EXTEND2_QP_UNROLL={toggle}",
                    "-include", str(ROOT / "include/swbwa_config.h"),
                    "-I", str(ROOT), "-I", str(ROOT / "include"),
                    "-I", str(ROOT / "tests/stubs"), str(source),
                    "-o", str(executable),
                ]
                print(f"Compiling actual ksw_extend2: {optimization} toggle={toggle}",
                      flush=True)
                subprocess.run(command, check=True)
                result = subprocess.run([str(executable)], stdout=subprocess.PIPE,
                                        check=True).stdout
                if canonical is None:
                    canonical = result
                if result != canonical:
                    raise AssertionError(f"extend2 output differs: {optimization}, {toggle}")
                print(f"PASS sha256={hashlib.sha256(result).hexdigest()}", flush=True)


if __name__ == "__main__":
    main()
