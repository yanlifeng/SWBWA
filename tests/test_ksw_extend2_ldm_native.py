"""Differential host tests of actual ksw_extend2 LDM placement and fallback.

Compile one snapshot of src/slave/ksw.c with explicit feature flags. The
reference uses heap storage and row-major QP construction (LDM=0, QP=0);
all builds disable SWBWA_ENABLE_CPE_KERNEL_OPT and enable only the requested
extend2 options. No Git history is required. The harness includes the actual
source unchanged and uses the real tracked
allocator, host LDM SDK stubs, and system/pool heap backing. It checks exact
outputs, EH/QP state, capacity boundaries, refusal fallback and ownership.
No source is edited. This is not a target performance test.
Run from any directory: python3 <repository>/tests/test_ksw_extend2_ldm_native.py
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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "malloc_wrap.h"

static void *test_ldm_alloc(unsigned long bytes, int site);
static void test_ldm_release(void *ptr, unsigned long bytes);
static void *test_realloc(void *ptr, size_t bytes)
{
#if SWBWA_TEST_POOL
    return wrap_realloc(ptr, bytes, __FILE__, __LINE__, __func__);
#else
    return realloc(ptr, bytes);
#endif
}
static void test_free(void *ptr)
{
#if SWBWA_TEST_POOL
    wrap_free(ptr, __FILE__, __LINE__, __func__);
#else
    free(ptr);
#endif
}
#define realloc test_realloc
#define free test_free
#define swbwa_ldm_alloc test_ldm_alloc
#define swbwa_ldm_release test_ldm_release
#include SWBWA_TEST_KSW_SOURCE
#undef realloc
#undef free
#undef swbwa_ldm_alloc
#undef swbwa_ldm_release

#ifndef SWBWA_KSW_EXTEND2_LDM_MAX_BYTES
#define SWBWA_KSW_EXTEND2_LDM_MAX_BYTES 0
#endif

_Static_assert(sizeof(eh_t) == 8, "EH footprint changed");
_Static_assert(sizeof(*swbwa_extend2_scratch.qp) == 1,
               "This LDM ablation tests byte QP storage");

int swbwa_test_cpe_id;
static int refuse_sdk;
static void *active_ldm;
static size_t active_bytes;
static unsigned long attempts, successes, releases, fallbacks;
static int current_qlen, current_m;
static const uint8_t *current_query;
static const int8_t *current_mat;
static uint64_t last_state_hash;

void *ldm_malloc(unsigned long bytes)
{
    void *ptr = refuse_sdk ? NULL : malloc(bytes);
    if (ptr) memset(ptr, 0xa5, bytes);
    return ptr;
}
void ldm_free(void *ptr, unsigned long bytes) { (void)bytes; free(ptr); }
void swbwa_cpe_fail(int code, long a, long b, long c)
{
    fprintf(stderr, "unexpected CPE error: %d %ld %ld %ld\n", code, a, b, c);
    abort();
}

static uint64_t check_state(const eh_t *eh, const int8_t *qp)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    const uint8_t *bytes = (const uint8_t *)eh;
    for (size_t i = 0; i < ((size_t)current_qlen + 1) * sizeof(eh_t); ++i)
        hash = (hash ^ bytes[i]) * UINT64_C(1099511628211);
    for (int row = 0; row < current_m; ++row) {
        for (int j = 0; j < current_qlen; ++j) {
            int8_t actual = qp[(size_t)row * current_qlen + j];
            assert(actual == current_mat[row * current_m + current_query[j]]);
            hash = (hash ^ (uint8_t)actual) * UINT64_C(1099511628211);
        }
    }
    return hash;
}

static void *test_ldm_alloc(unsigned long bytes, int site)
{
    void *ptr;
    assert(site == 8 && !active_ldm);
    assert(bytes == ((size_t)current_qlen + 1) * sizeof(eh_t) +
                    (size_t)current_qlen * current_m);
    assert(bytes <= SWBWA_KSW_EXTEND2_LDM_MAX_BYTES);
    ++attempts;
    ptr = swbwa_ldm_alloc(bytes, site);
    if (ptr) {
        ++successes;
        active_ldm = ptr;
        active_bytes = bytes;
        assert((uintptr_t)ptr % _Alignof(eh_t) == 0);
    }
    return ptr;
}

static void test_ldm_release(void *ptr, unsigned long bytes)
{
    assert(ptr == active_ldm && bytes == active_bytes);
    last_state_hash = check_state(ptr, (int8_t *)((eh_t *)ptr + current_qlen + 1));
    memset(ptr, 0xdd, bytes);
    active_ldm = NULL;
    ++releases;
    swbwa_ldm_release(ptr, bytes);
}

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
    uint8_t saved_query[800], saved_target[800];
    int8_t saved_mat[400];
    size_t live_qp = (size_t)qlen * m;
    size_t live_eh = ((size_t)qlen + 1) * sizeof(eh_t);
    size_t live_bytes = live_eh + live_qp;
    unsigned mode = cases % 6;
    unsigned long attempts_before = attempts, successes_before = successes;
    unsigned long releases_before = releases;
    swbwa_extend2_scratch_t saved_scratch;
    size_t held_bytes = mode == 1 ? 17520 : mode == 2 ? 40960 - live_bytes :
                        mode == 3 ? 40960 - live_bytes + 1 : mode == 5 ? 40960 : 0;
    void *held = NULL;
    int eligible = live_bytes <= SWBWA_KSW_EXTEND2_LDM_MAX_BYTES;
    int should_succeed = eligible && mode != 3 && mode != 4 && mode != 5;

    assert(qlen > 0 && qlen <= 800 && tlen >= 0 && tlen <= 800);
    assert(m > 0 && m <= 20 && ed > 0 && ei > 0);
    assert(live_bytes <= 40960 && swbwa_ldm_outstanding() == 0);
    refuse_sdk = 0;
    if (held_bytes) {
        held = swbwa_ldm_alloc(held_bytes, 4);
        assert(held);
    }
    swbwa_ldm_begin_batch();
    refuse_sdk = mode == 4;
    current_qlen = qlen; current_m = m;
    current_query = query; current_mat = mat;
    memcpy(saved_query, query, qlen);
    memcpy(saved_target, target, tlen);
    memcpy(saved_mat, mat, (size_t)m * m);
    if (swbwa_extend2_scratch.qp)
        memset(swbwa_extend2_scratch.qp, 0xa5, swbwa_extend2_scratch.qp_capacity);
    if (swbwa_extend2_scratch.eh)
        memset(swbwa_extend2_scratch.eh, 0x5a, swbwa_extend2_scratch.eh_capacity);
    saved_scratch = swbwa_extend2_scratch;
    for (int i = 0; i < 5; ++i) ptr[i] = null_mask & (1u << i) ? NULL : &out[i];
    score = ksw_extend2(qlen, query, tlen, target, m, mat,
                       od, ed, oi, ei, w, bonus, zdrop, h0,
                       ptr[0], ptr[1], ptr[2], ptr[3], ptr[4]);
    ++cases;
    assert(attempts - attempts_before == (unsigned)eligible);
    assert(successes - successes_before == (unsigned)should_succeed);
    assert(releases - releases_before == (unsigned)should_succeed);
    assert(!active_ldm && swbwa_ldm_outstanding() == (long)held_bytes);
    assert(swbwa_ldm_peak() <= 40960);
    assert(swbwa_ldm_refusals() == (eligible && !should_succeed));
    if (should_succeed) {
        assert(swbwa_extend2_scratch.eh == saved_scratch.eh);
        assert(swbwa_extend2_scratch.qp == saved_scratch.qp);
        assert(swbwa_extend2_scratch.eh_capacity == saved_scratch.eh_capacity);
        assert(swbwa_extend2_scratch.qp_capacity == saved_scratch.qp_capacity);
    } else {
        ++fallbacks;
        last_state_hash = check_state(swbwa_extend2_scratch.eh,
                                      swbwa_extend2_scratch.qp);
    }
    qp_bytes += live_qp;
    for (size_t i = should_succeed ? 0 : live_qp;
         i < swbwa_extend2_scratch.qp_capacity; ++i)
        assert((uint8_t)swbwa_extend2_scratch.qp[i] == 0xa5);
    for (size_t i = should_succeed ? 0 : live_eh;
         i < swbwa_extend2_scratch.eh_capacity; ++i)
        assert(((uint8_t *)swbwa_extend2_scratch.eh)[i] == 0x5a);
    for (int i = 0; i < 5; ++i)
        if (null_mask & (1u << i)) assert(out[i] == -777);
    assert(memcmp(query, saved_query, qlen) == 0);
    assert(memcmp(target, saved_target, tlen) == 0);
    assert(memcmp(mat, saved_mat, (size_t)m * m) == 0);
    printf("%lu %d %d %d %d %d %d %016llx\n", cases, score,
           out[0], out[1], out[2], out[3], out[4],
           (unsigned long long)last_state_hash);
    refuse_sdk = 0;
    if (held) swbwa_ldm_release(held, held_bytes);
    assert(swbwa_ldm_outstanding() == 0);
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
                           63, 64, 75, 127, 128, 129, 149, 150, 151, 255, 256,
                           313, 314, 315, 511, 628, 629, 630, 799};
    const int alphabets[] = {1, 2, 4, 5, 6, 20};
    uint8_t query[800], target[800];
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
#if SWBWA_TEST_POOL
    char *pool = malloc(16 << 20);
    assert(pool);
    set_big_buffer(pool, 16 << 20);
#endif
    test_short_exhaustive();
    test_lengths_matrices_and_fallback();
    fprintf(stderr, "LDM_CAP=%d QP_UNROLL=%d POOL=%d: %lu actual extend2 calls, "
            "%lu QP bytes checked; attempts=%lu successes=%lu releases=%lu heap=%lu\n",
            SWBWA_KSW_EXTEND2_LDM_MAX_BYTES, SWBWA_KSW_EXTEND2_QP_UNROLL,
            SWBWA_TEST_POOL, cases, qp_bytes, attempts, successes, releases, fallbacks);
    assert(successes == releases && swbwa_ldm_outstanding() == 0);
    test_free(swbwa_extend2_scratch.qp);
    test_free(swbwa_extend2_scratch.eh);
#if SWBWA_TEST_POOL
    free(pool);
#endif
    return 0;
}
'''


def main():
    compiler = shlex.split(os.environ.get("CC", "gcc"))
    contents = (ROOT / "src/slave/ksw.c").read_bytes()
    with tempfile.TemporaryDirectory(prefix="swbwa-extend2-ldm-") as temporary:
        path = Path(temporary)
        source = path / "test.c"
        source.write_text(HARNESS, encoding="ascii")
        snapshot = path / "ksw_snapshot.c"
        snapshot.write_bytes(contents)
        print(f"shared source sha256={hashlib.sha256(contents).hexdigest()}",
              flush=True)
        canonical = None
        matrix = [(optimization, pool, qp, cap)
                  for optimization in ("-O2", "-O3") for pool in (0, 1)
                  for cap in (0, 4096, 8192) for qp in (0, 1)]
        for optimization, pool, toggle, cap in matrix:
            name = "reference" if cap == 0 and toggle == 0 else "candidate"
            executable = path / f"test-{optimization[1:]}-{pool}-{toggle}-{name}-{cap}"
            command = compiler + [
                "-std=gnu11", optimization, "-g", "-Wall", "-Wextra",
                "-Wno-unused-parameter", "-Wno-sign-compare",
                "-Werror=implicit-function-declaration",
                "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections",
                "-fsanitize=address,undefined", "-fno-sanitize-recover=all",
                "-fno-omit-frame-pointer", "-fno-pie", "-no-pie",
                "-DSWBWA_ENABLE_CPE_KERNEL_OPT=0",
                "-DSWBWA_ENABLE_CPE_PROFILE=0",
                "-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=0",
                "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                "-DSWBWA_ENABLE_MATESW_DUAL_FORWARD=0",
                f"-DSWBWA_TEST_POOL={pool}",
                f"-DSWBWA_KSW_EXTEND2_LDM_MAX_BYTES={cap}",
                f'-DSWBWA_TEST_KSW_SOURCE="{snapshot}"',
                f"-DSWBWA_KSW_EXTEND2_QP_UNROLL={toggle}",
                "-include", str(ROOT / "include/swbwa_config.h"),
                "-I", str(ROOT), "-I", str(ROOT / "src/slave"),
                "-I", str(ROOT / "include"),
                "-I", str(ROOT / "tests/stubs"), str(source),
                str(ROOT / "src/slave/malloc_wrap.c"), "-lm", "-o", str(executable),
            ]
            print(f"Compiling {name}: {optimization} pool={pool} qp={toggle} cap={cap}",
                  flush=True)
            subprocess.run(command, check=True)
            result = subprocess.run([str(executable)], stdout=subprocess.PIPE,
                                    check=True).stdout
            if canonical is None:
                assert name == "reference", "Run the heap/QP-off reference first"
                canonical = result
            if result != canonical:
                raise AssertionError(f"extend2 output differs: {name}, {optimization}, "
                                     f"pool={pool}, qp={toggle}, cap={cap}")
            print(f"PASS sha256={hashlib.sha256(result).hexdigest()}", flush=True)
        print(f"PASS: all {len(matrix)} variants match the same-source "
              "heap/QP-off reference (SWBWA_ENABLE_CPE_KERNEL_OPT=0)", flush=True)


if __name__ == "__main__":
    main()
