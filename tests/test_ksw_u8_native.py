"""Compare actual whole int32_16 ksw_u8 with fused/XOR independently off/on.

Run: python3 tests/test_ksw_u8_native.py. Uses the existing GCC-vector shim,
unchanged extracted production functions, ASan and UBSan. This is not a
Sunway compile, performance test, alignment-wrapper test, or SAM fingerprint.
"""

import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from test_ksw_i16_native import PREAMBLE, ROOT, SIMD, extract


PROFILE = r'''
static int rows, stripe;
static uint64_t lazy;
static void record_matesw_ksw_work(int r, int s, uint64_t l) {
    rows = r; stripe = s; lazy = l;
}
static_assert(SWBWA_KSW_U8_MODE == SWBWA_KSW_U8_INT32_16,
              "this test requires the int32_16 u8 kernel");
'''

HARNESS = r'''
static uint32_t seed = 7963242;
static unsigned cases, skipped_sums, saturated, suboptimal, stopped;
static uint32_t random_word() {
    seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
    return seed;
}
static uint64_t hash_bytes(const void *ptr, size_t bytes) {
    const unsigned char *p = (const unsigned char *)ptr;
    uint64_t h = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < bytes; ++i)
        h = (h ^ p[i]) * UINT64_C(1099511628211);
    return h;
}
static void check_profile(const kswq_t *q, const std::vector<uint8_t> &query,
                          int m, const std::vector<int8_t> &mat) {
    int minimum = 127, maximum = 0;
    for (int v : mat) {
        minimum = std::min(minimum, v); maximum = std::max(maximum, v);
    }
    assert(q->shift == (uint8_t)(-minimum));
    assert(q->max == maximum && q->mdiff == (uint8_t)(maximum + q->shift));
    for (int a = 0; a < m; ++a) {
        for (int i = 0; i < q->slen; ++i) {
            int words[16]; simd_store(q->qp[a * q->slen + i].words, words);
            for (int lane = 0; lane < 16; ++lane) {
                int pos = i + lane * q->slen;
                int expected = (pos < q->qlen ? mat[a * m + query[pos]] : 0)
                    + q->shift;
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
        ++skipped_sums; return;
    }
    std::vector<uint8_t> query(qlen), target(tlen);
    std::vector<int8_t> mat(m * m);
    for (int a = 0; a < m; ++a)
        for (int b = 0; b < m; ++b)
            mat[a * m + b] = scoring == 0 ? (a == b ? 1 : -4)
                : scoring == 1 ? (a == b ? 2 : -8)
                : scoring == 2 ? (a == b ? 1 : 0)
                : scoring == 4 ? 0 : (int)(random_word() % 256) - 128;
    // A positive maximum avoids the existing XSUBO radius division by zero.
    if (scoring == 3) mat[0] = 127;
    for (int i = 0; i < qlen; ++i)
        query[i] = pattern == 0 || pattern == 5 ? 0 : pattern == 1 ? i % m
            : random_word() % m;
    for (int i = 0; i < tlen; ++i)
        target[i] = pattern == 0 ? 0 : pattern == 1 ? i % m
            : pattern == 2 ? query[i % qlen]
            : pattern == 5 ? (i < 30 || (i >= 120 && i < 150) ? 0 : 1)
            : random_word() % m;
    const std::vector<uint8_t> saved_query = query;
    const std::vector<int8_t> saved_mat = mat;
    kswq_t *q = ksw_qinit_impl(1, qlen, query.data(), m, mat.data(), 0);
    assert(q->size == 1 && q->slen == (qlen + 15) / 16 && !q->in_ldm);
    check_profile(q, query, m, mat);
    // Define H1 even for zero target rows, without changing kernel inputs H0/E.
    std::memset(q->H0, 0, 4 * q->slen * sizeof(__m128i));
    for (int repeat = 0; repeat < (reuse ? 2 : 1); ++repeat) {
        const std::vector<uint8_t> saved_target = target;
        rows = stripe = -1; lazy = UINT64_MAX;
        kswr_t r = ksw_u8(q, tlen, target.data(), od, ed, oi, ei, xtra);
        assert(query == saved_query && target == saved_target && mat == saved_mat);
        check_profile(q, query, m, mat);
        assert(rows >= 0 && rows <= tlen && stripe == q->slen);
        assert(lazy <= (uint64_t)rows * q->slen * 16);
        if (r.score == 255) {
            ++saturated;
            assert(r.qe == -1 && r.score2 == -1 && r.te2 == -1);
        }
        if (r.score2 >= 0) ++suboptimal;
        if (rows < tlen) ++stopped;
        assert(r.tb == -1 && r.qb == -1);
        if (pattern == 0 && scoring == 0 && m == 5 && od == 6 && ed == 1 &&
            oi == 6 && ei == 1 && xtra == 0) {
            int length = std::min(qlen, tlen);
            assert(r.score == (length >= 251 ? 255 : length));
            assert(r.te == std::min(length, 251) - 1);
            assert(r.qe == (length >= 251 ? -1 : tlen ? length - 1 : 0));
        }
        if (pattern == 5) {
            assert(r.score == 30 && r.te == 29 && r.qe == 29);
            assert(r.score2 == 30 && r.te2 == 149);
        }
        std::printf("%u %d %d %d %d %d %d %d %d %d %llu %016llx %016llx\n",
                    ++cases, r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb,
                    rows, stripe, (unsigned long long)lazy,
                    (unsigned long long)hash_bytes(q->qp,
                                                  m * q->slen * sizeof(__m128i)),
                    (unsigned long long)hash_bytes(q->H0,
                                                  4 * q->slen * sizeof(__m128i)));
        std::reverse(target.begin(), target.end());
    }
    std::free(q);
}
int main() {
    const int lengths[] = {1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65,
                           75, 127, 128, 129, 149, 150, 151, 250, 251, 255,
                           256, 257, 300, 511};
    const int flags[] = {0, KSW_XSUBO, KSW_XSUBO | 1, KSW_XSUBO | 20,
                         KSW_XSUBO | 255, KSW_XSTOP | 1, KSW_XSTOP | 249,
                         KSW_XSTOP | KSW_XSUBO | 30,
                         KSW_XSTOP | KSW_XSUBO | 65535};
    for (int qlen : lengths)
        for (int tlen : {0, 1, qlen, qlen + 17})
            for (int pattern = 0; pattern < 4; ++pattern)
                for (int xtra : flags)
                    run_case(qlen, tlen, 5, pattern, pattern,
                             6, 1, 6, 1, xtra, true);
    run_case(30, 220, 5, 5, 0, 6, 1, 6, 1, KSW_XSUBO | 20);
    run_case(33, 47, 5, 0, 4, 6, 1, 6, 1, 0);
    for (int n = 0; n < 2000; ++n) {
        int args[10];
        for (int &v : args) v = random_word() & INT_MAX;
        run_case(1 + args[0] % 300, args[1] % 400, 1 + args[2] % 8,
                 args[3] % 4, args[4] % 4, args[5] % 32, args[6] % 5,
                 args[7] % 32, args[8] % 5, flags[args[9] % 9]);
    }
    // All 256x256 raw opening/extension pairs, including signed-byte narrowing.
    for (int opening = 0; opening < 256; ++opening)
        for (int extension = 0; extension < 256; ++extension)
            run_case(1 + opening % 33, 1 + extension % 9, 5,
                     opening % 4, extension % 4, opening, extension,
                     255 - extension, 255 - opening, flags[extension % 9]);
    const int penalties[] = {INT_MIN, INT_MIN + 1, -65537, -65536, -32769,
        -32768, -257, -256, -129, -128, -1, 0, 1, 127, 128, 255, 256,
        32767, 32768, 65535, 65536, 65537, INT_MAX - 1, INT_MAX};
    for (int opening : penalties)
        for (int extension : penalties)
            run_case(17, 9, 5, 2, 3, opening, extension, 1, 2, KSW_XSUBO | 1);
    assert(saturated && suboptimal && stopped);
    std::fprintf(stderr, "%u whole-kernel calls; saturated=%u suboptimal=%u "
                 "early-exit=%u; %u undefined raw sums excluded\n",
                 cases, saturated, suboptimal, stopped, skipped_sums);
}
'''


def main():
    source = (ROOT / "src/slave/ksw.c").read_text()
    extracted = extract(source, "const kswr_t g_defr =", "\n#ifndef SWBWA_KSW_LDM")
    extracted += extract(source, "static kswq_t *ksw_qinit_impl(",
                         "\nkswq_t *ksw_qinit(")
    extracted += extract(source, "SWBWA_MATESW_HOT\nkswr_t ksw_u8(",
                         "\nSWBWA_MATESW_HOT\nkswr_t ksw_i16(")
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    environment = dict(os.environ, ASAN_OPTIONS="detect_leaks=1")
    with tempfile.TemporaryDirectory(prefix="swbwa-ksw-u8-") as temporary:
        path = Path(temporary)
        (path / "simd.h").write_text(SIMD, encoding="ascii")
        (path / "test.cpp").write_text(
            PREAMBLE + PROFILE + extracted + '\n#line 1 "test-harness.cpp"\n' + HARNESS,
            encoding="ascii")
        baseline = None
        for fused, xor in ((0, 0), (1, 0), (0, 1), (1, 1)):
            name = f"fused{fused}-xor{xor}"
            executable = path / name
            print(f"Actual ksw_u8 + scalar_sse.h: {name}", flush=True)
            subprocess.run(compiler + [
                "-std=c++11", "-O2", "-g0", "-Wall", "-Wextra", "-Werror",
                "-Wno-psabi", "-fsanitize=address,undefined", "-no-pie",
                "-fno-sanitize-recover=all", "-DSWBWA_ENABLE_CPE_PROFILE=1",
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
            print(f"  {len(output.splitlines())} results/QP/DP/profile records match; "
                  f"SHA256={hashlib.sha256(output).hexdigest()}", flush=True)
    print("PASS: baseline, fused-only, XOR-only, fused+XOR; ASan+UBSan")


if __name__ == "__main__":
    main()
