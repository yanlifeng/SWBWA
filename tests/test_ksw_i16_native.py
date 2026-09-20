"""Differentially compile actual ksw_i16 and scalar_sse.h with a SIMD shim.

Run: python3 tests/test_ksw_i16_native.py (CXX defaults to g++).
The initializer and whole kernel are extracted unchanged, not reimplemented.
This tests CPE integer dataflow, not Sunway lowering or full SAM fingerprints.
Raw int gap-open + extension overflow is baseline UB and deliberately excluded.
Generated sources, binaries, and comparison output stay in a temporary directory.
"""

import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from test_ksw_fused_gap_native import SIMD


ROOT = Path(__file__).resolve().parents[1]

PREAMBLE = r'''
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "swbwa_config.h"
#include "src/slave/scalar_sse.h"
#include "src/slave/ksw.h"

// Select the actual CPE branches without changing the native compiler ABI.
#undef __SSE2__
#undef __ARM_NEON
#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#define SWBWA_MATESW_HOT __attribute__((optimize("O3")))
#define SWBWA_KSW_I16_SHIFT_LEFT_ONE(v) swbwa_i16_shift_left_lane(v)
#define SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES (16 << 10)
enum { SWBWA_KSW_U8_LANES = 16 };
void *swbwa_ldm_alloc(unsigned long, int) { std::abort(); }
static_assert(SWBWA_KSW_I16_MODE == SWBWA_KSW_I16_INT32_8,
              "this test exercises the int32_8 CPE i16 kernel");
'''

HARNESS = r'''
static uint32_t seed = 7963242;
static unsigned cases, skipped_sums;
static uint32_t random_word() {
    seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
    return seed;
}
static uint64_t state_hash(const void *ptr, size_t size) {
    const unsigned char *bytes = (const unsigned char *)ptr;
    uint64_t h = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < size; ++i)
        h = (h ^ bytes[i]) * UINT64_C(1099511628211);
    return h;
}
static void check_profile(kswq_t *q, const std::vector<uint8_t> &query,
                          int m, const std::vector<int8_t> &mat) {
    for (int a = 0; a < m; ++a) {
        for (int i = 0; i < q->slen; ++i) {
            int words[16]; simd_store(q->qp[a * q->slen + i].words, words);
            for (int lane = 0; lane < 16; ++lane) {
                int pos = i + lane * q->slen;
                int expected = lane < 8 && pos < q->qlen
                    ? mat[a * m + query[pos]] : 0;
                assert(words[lane] == expected);
            }
        }
    }
}
static void run_case(int qlen, int tlen, int m, int pattern, int scoring,
                     int od, int ed, int oi, int ei, int xtra,
                     bool reuse = false) {
    int64_t sd = (int64_t)od + ed, si = (int64_t)oi + ei;
    if (sd < INT_MIN || sd > INT_MAX || si < INT_MIN || si > INT_MAX) {
        ++skipped_sums;
        return;
    }
    std::vector<uint8_t> query(qlen), target(tlen);
    std::vector<int8_t> mat(m * m);
    for (int a = 0; a < m; ++a)
        for (int b = 0; b < m; ++b)
            mat[a * m + b] = scoring == 0 ? (a == b ? 2 : -8)
                : scoring == 1 ? (a == b ? 127 : -128)
                : scoring == 2 ? (a == b ? 1 : 0)
                : (int)(random_word() % 256) - 128;
    // XSUBO's existing radius calculation requires a positive matrix maximum.
    mat[0] = scoring == 1 ? 127 : 2;
    for (int i = 0; i < qlen; ++i)
        query[i] = pattern == 0 ? 0 : pattern == 1 ? i % m
            : random_word() % m;
    for (int i = 0; i < tlen; ++i)
        target[i] = pattern == 0 ? 0 : pattern == 1 ? i % m
            : pattern == 2 ? query[i % qlen] : random_word() % m;
    const std::vector<uint8_t> saved_query = query;
    const std::vector<int8_t> saved_mat = mat;
    kswq_t *q = ksw_qinit_impl(2, qlen, query.data(), m, mat.data(), 0);
    assert(q->slen == (qlen + 7) / 8 && !q->in_ldm);
    check_profile(q, query, m, mat);
    // H1 is not initialized by the kernel when tlen==0. Initialize all scratch
    // here so its final-state checksum is deterministic without changing H0/E.
    std::memset(q->H0, 0, 4 * q->slen * sizeof(__m128i));
    for (int repeat = 0; repeat < (reuse ? 2 : 1); ++repeat) {
        const std::vector<uint8_t> saved_target = target;
        kswr_t r = ksw_i16(q, tlen, target.data(), od, ed, oi, ei, xtra);
        assert(query == saved_query && mat == saved_mat && target == saved_target);
        check_profile(q, query, m, mat);
        for (int i = 0; i < 4 * q->slen; ++i) {
            int words[16]; simd_store(q->H0[i].words, words);
            for (int lane = 8; lane < 16; ++lane) assert(words[lane] == 0);
        }
        if (pattern == 0 && scoring == 0 && od == 12 && ed == 2 &&
            oi == 12 && ei == 2 && xtra == 0) {
            assert(r.score == 2 * std::min(qlen, tlen));
            assert(r.te == std::min(qlen, tlen) - 1);
            assert(r.qe == (tlen ? std::min(qlen, tlen) - 1 : 0));
        }
        if (pattern == 0 && scoring == 1 && qlen == 300 && tlen == 330 &&
            od == 12 && ed == 2 && oi == 12 && ei == 2 && xtra == 0)
            assert(r.score == 32767);
        // Compare every kswr_t field directly; also fingerprint all DP words.
        std::printf("%u %d %d %d %d %d %d %d %016llx\n", ++cases,
                    r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb,
                    (unsigned long long)state_hash(q->H0,
                                                  4 * q->slen * sizeof(__m128i)));
        std::reverse(target.begin(), target.end());
    }
    std::free(q);
}
int main() {
    const int lengths[] = {1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33,
                           75, 127, 128, 129, 149, 150, 151, 255, 256, 300, 511};
    const int flags[] = {0, KSW_XSUBO, KSW_XSUBO | 1, KSW_XSUBO | 255,
                         KSW_XSTOP | 1, KSW_XSTOP | 255,
                         KSW_XSTOP | KSW_XSUBO | 30,
                         KSW_XSTOP | KSW_XSUBO | 65535};
    for (int qlen : lengths) {
        for (int tlen : {0, 1, qlen, qlen + 17}) {
            for (int pattern = 0; pattern < 4; ++pattern) {
                for (int xtra : flags)
                    run_case(qlen, tlen, 5, pattern, pattern,
                             12, 2, 12, 2, xtra, true);
            }
        }
    }
    run_case(300, 330, 5, 0, 1, 12, 2, 12, 2, 0);
    for (int n = 0; n < 2000; ++n)
        run_case(1 + random_word() % 180, random_word() % 240,
                 1 + random_word() % 8, random_word() % 4, random_word() % 4,
                 random_word() % 32, random_word() % 5,
                 random_word() % 32, random_word() % 5,
                 flags[random_word() % 8]);
    // Every signed-16 bit pattern, with distinct insertion/deletion penalties.
    for (int raw = 0; raw <= 65535; ++raw)
        run_case(1 + raw % 17, 1 + raw % 7, 5, raw % 4, raw % 4,
                 1, raw, 17, 65535 - raw, flags[raw % 8]);
    const int penalties[] = {INT_MIN, INT_MIN + 1, -65537, -65536, -32769,
        -32768, -1, 0, 1, 127, 128, 255, 256, 32767, 32768, 65535, 65536,
        65537, INT_MAX - 1, INT_MAX};
    for (int opening : penalties)
        for (int extension : penalties)
            run_case(17, 9, 5, 2, 1, opening, extension, 1, 2,
                     KSW_XSUBO | 1);
    std::fprintf(stderr, "%u whole-kernel calls; %u undefined raw sums excluded\n",
                 cases, skipped_sums);
}
'''


def extract(source, start, end):
    """Fail closed on moved markers; emit verbatim code with original line IDs."""
    assert source.count(start) == 1, start
    begin = source.index(start)
    finish = source.index(end, begin + len(start))
    return (f'#line {source.count(chr(10), 0, begin) + 1} "src/slave/ksw.c"\n'
            + source[begin:finish] + "\n")


def main():
    source = (ROOT / "src/slave/ksw.c").read_text()
    extracted = extract(source, "const kswr_t g_defr =", "\n#ifndef SWBWA_KSW_LDM")
    extracted += extract(source, "static kswq_t *ksw_qinit_impl(",
                         "\nkswq_t *ksw_qinit(")
    extracted += extract(source, "SWBWA_MATESW_HOT\nkswr_t ksw_i16(",
                         "\n#if SWBWA_ENABLE_MATESW_DUAL_FORWARD &&")
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    environment = dict(os.environ, ASAN_OPTIONS="detect_leaks=1")
    with tempfile.TemporaryDirectory(prefix="swbwa-ksw-i16-") as temporary:
        path = Path(temporary)
        (path / "simd.h").write_text(SIMD, encoding="ascii")
        (path / "test.cpp").write_text(
            PREAMBLE + extracted + '\n#line 1 "test-harness.cpp"\n' + HARNESS,
            encoding="ascii")
        baseline = None
        for fused, xor in ((0, 0), (1, 0), (0, 1), (1, 1)):
            name = f"fused{fused}-xor{xor}"
            executable = path / name
            print(f"Actual ksw_i16 + scalar_sse.h: {name}", flush=True)
            subprocess.run(compiler + [
                "-std=c++11", "-O2", "-g0", "-Wall", "-Wextra", "-Werror",
                "-Wno-psabi", "-fsanitize=address,undefined", "-no-pie",
                "-fno-sanitize-recover=all", "-DSWBWA_ENABLE_CPE_PROFILE=0",
                "-DSWBWA_KSW_I16_MODE=2", "-DSWBWA_KSW_U8_MODE=1",
                f"-DSWBWA_KSW_FUSED_GAP_UPDATE={fused}",
                f"-DSWBWA_KSW_XOR_SELECT={xor}",
                "-I", str(path), "-I", str(ROOT / "include"), "-I", str(ROOT),
                str(path / "test.cpp"), "-o", str(executable),
            ], check=True)
            output = subprocess.check_output([str(executable)], env=environment)
            if baseline is None:
                baseline = output
            elif output != baseline:
                for i, (old, new) in enumerate(zip(baseline.splitlines(),
                                                 output.splitlines()), 1):
                    if old != new:
                        raise AssertionError(f"{name} case {i}: {old!r} != {new!r}")
                raise AssertionError(f"{name}: output length differs")
            print(f"  {len(output.splitlines())} results/DP states match; "
                  f"SHA256={hashlib.sha256(output).hexdigest()}", flush=True)
    print("PASS: baseline, fused-only, XOR-only, fused+XOR; ASan+UBSan")


if __name__ == "__main__":
    main()
