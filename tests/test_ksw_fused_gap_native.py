"""Compile the unmodified scalar_sse.h with a GCC-vector Sunway SIMD shim.

Run with python3 tests/test_ksw_fused_gap_native.py; CXX defaults to g++.
The shim only implements integer word operations. It tests actual header
functions, not Sunway instruction lowering, performance, or SAM fingerprints.
All generated sources and binaries live in an automatically removed temp dir.
"""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]

SIMD = r'''
#ifndef SWBWA_TEST_GAP_SIMD_H
#define SWBWA_TEST_GAP_SIMD_H
#include <cassert>
#include <cstdint>
#include <cstring>
#include <initializer_list>

typedef int32_t word16 __attribute__((vector_size(64)));
typedef uint32_t uword16 __attribute__((vector_size(64)));

/* GCC does not accept Sunway's implicit scalar-to-vector broadcasts. Only
 * that syntax is adapted here; primitive arithmetic uses real GCC vectors.
 */
struct alignas(64) intv16 {
    word16 lanes;
    intv16() = default;
    intv16(int32_t value) {
        for (int i = 0; i < 16; ++i) lanes[i] = value;
    }
    intv16(std::initializer_list<int32_t> values) {
        assert(values.size() == 16);
        int i = 0;
        for (int32_t value : values) lanes[i++] = value;
    }
};
static_assert(sizeof(intv16) == 64, "wrong vector size");

static inline intv16 simd_vaddw(intv16 a, intv16 b) {
    intv16 r;
    r.lanes = (word16)((uword16)a.lanes + (uword16)b.lanes);
    return r;
}
static inline intv16 simd_vsubw(intv16 a, intv16 b) {
    intv16 r;
    r.lanes = (word16)((uword16)a.lanes - (uword16)b.lanes);
    return r;
}
static inline intv16 simd_vandw(intv16 a, intv16 b) {
    intv16 r; r.lanes = a.lanes & b.lanes; return r;
}
static inline intv16 simd_vxorw(intv16 a, intv16 b) {
    intv16 r; r.lanes = a.lanes ^ b.lanes; return r;
}
static inline intv16 simd_vbisw(intv16 a, intv16 b) {
    intv16 r; r.lanes = a.lanes | b.lanes; return r;
}
static inline intv16 simd_vcmpltw(intv16 a, intv16 b) {
    intv16 r, one = 1;
    /* Sunway comparison lanes are 0/1, unlike GCC's 0/-1. */
    r.lanes = (a.lanes < b.lanes) & one.lanes;
    return r;
}
static inline void simd_load(intv16 &r, const int *p) {
    std::memcpy(&r.lanes, p, 64);
}
static inline void simd_store(intv16 a, int *p) {
    std::memcpy(p, &a.lanes, 64);
}
static inline intv16 simd_sllx(intv16 a, int bits) {
    assert(bits >= 0 && bits <= 512 && bits % 32 == 0);
    intv16 r = 0;
    for (int i = bits / 32; i < 16; ++i)
        r.lanes[i] = a.lanes[i - bits / 32];
    return r;
}
static inline int simd_vextw15(intv16 a) { return a.lanes[15]; }
#endif
'''

HARNESS = r'''
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "swbwa_config.h"
#include "src/slave/scalar_sse.h"

#if !SWBWA_KSW_FUSED_GAP_UPDATE || SWBWA_ENABLE_FLOAT16_VECTOR || SWBWA_ENABLE_PACKED_INT8
#error "this test requires fused integer word helpers"
#endif

static uint64_t cases, lane_checks, proposal_cases, skipped_undefined_sums;
static uint32_t rng_state = 7963242;

static uint32_t random_word() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}
static int32_t signed_word(uint32_t bits) {
    int32_t value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}
static int32_t subtract_word(int32_t a, int32_t b) {
    return signed_word((uint32_t)a - (uint32_t)b);
}
static int32_t narrow_signed(int raw, int bits) {
    uint32_t mask = (UINT32_C(1) << bits) - 1;
    int32_t low = (uint32_t)raw & mask;
    return low >= (1 << (bits - 1)) ? low - (1 << bits) : low;
}
static __m128i load(const int *values) {
    __m128i r;
    simd_load(r.words, values);
    return r;
}
static void equal(__m128i a, __m128i b, const char *phase) {
    int av[16], bv[16];
    simd_store(a.words, av);
    simd_store(b.words, bv);
    for (int i = 0; i < 16; ++i) {
        ++lane_checks;
        if (av[i] != bv[i]) {
            std::fprintf(stderr, "%s case=%llu lane=%d: %d != %d\n",
                         phase, (unsigned long long)cases, i, av[i], bv[i]);
            std::abort();
        }
    }
}
static __m128i check(bool i16, __m128i gap, __m128i ext,
                     __m128i h, __m128i open_ext) {
    /* These are the original expressions, using the actual header helpers. */
    __m128i original = i16
        ? _mm_max_epi16(_mm_subs_epu16(gap, ext), _mm_subs_epu16(h, open_ext))
        : _mm_max_epu8(_mm_subs_epu8(gap, ext), _mm_subs_epu8(h, open_ext));
    __m128i fused = i16
        ? swbwa_ksw_i16_gap_update(gap, ext, h, open_ext)
        : swbwa_ksw_gap_update_words(gap, ext, h, open_ext);
    int gv[16], ev[16], hv[16], ov[16], expected[16];
    simd_store(gap.words, gv); simd_store(ext.words, ev);
    simd_store(h.words, hv); simd_store(open_ext.words, ov);
    for (int i = 0; i < 16; ++i) {
        expected[i] = i16 && i >= 8 ? 0
            : std::max(0, std::max(subtract_word(gv[i], ev[i]),
                                   subtract_word(hv[i], ov[i])));
    }
    ++cases;
    equal(original, load(expected), "original/oracle");
    equal(fused, original, "fused/original");
    return fused;
}

static void test_shim() {
    intv16 lo = INT_MIN, hi = INT_MAX;
    assert(simd_vcmpltw(lo, hi).lanes[0] == 1);
    assert(simd_vcmpltw(hi, lo).lanes[0] == 0);
    assert(simd_vcmpltw(hi, hi).lanes[0] == 0);
    assert(simd_vsubw(lo, 1).lanes[0] == INT_MAX);
    assert(simd_vaddw(hi, 1).lanes[0] == INT_MIN);
    int values[16];
    for (int i = 0; i < 16; ++i) values[i] = i + 1;
    intv16 a; simd_load(a, values);
    a = simd_sllx(a, 32);
    assert(a.lanes[0] == 0 && simd_vextw15(a) == 15);
}

static void test_u8_differences() {
    for (int a = -255; a <= 255; ++a) {
        for (int b = -255; b <= 255; ++b) {
            __m128i g, e, h, o;
            g.words = std::max(a, 0); e.words = std::max(-a, 0);
            h.words = std::max(b, 0); o.words = std::max(-b, 0);
            check(false, g, e, h, o);
        }
    }
}

static void test_penalty_narrowing() {
    int state[16] = {INT16_MIN, -1, 0, 1, 255, 32767, 32768, 65535,
                    17, -3, 90000, -65535, 0, 1, 32767, -32768};
    __m128i g = load(state), h = load(state);
    for (int raw = 0; raw <= 65535; ++raw) {
        __m128i e = _mm_set1_epi16(raw);
        __m128i o = _mm_set1_epi16(65535 - raw);
        int expected[16];
        for (int i = 0; i < 16; ++i)
            expected[i] = i < 8 ? narrow_signed(raw, 16) : 0;
        equal(e, load(expected), "i16 broadcast truncation");
        check(true, g, e, h, o);
    }
    for (int raw = -512; raw <= 767; ++raw) {
        __m128i e = _mm_set1_epi8(raw);
        int expected[16];
        for (int i = 0; i < 16; ++i) expected[i] = narrow_signed(raw, 8);
        equal(e, load(expected), "u8 broadcast truncation");
        check(false, g, e, h, _mm_set1_epi8(raw + 17));
    }
    const int penalties[] = {INT_MIN, INT_MIN + 1, -131073, -65537, -65536,
        -32769, -32768, -1, 0, 1, 127, 128, 255, 256, 32766, 32767,
        32768, 65535, 65536, 65537, 131071, INT_MAX - 1, INT_MAX};
    for (int opening : penalties) {
        for (int extension : penalties) {
            int64_t sum = (int64_t)opening + extension;
            /* ksw_i16 adds int arguments before the int16_t conversion. */
            if (sum < INT_MIN || sum > INT_MAX) {
                ++skipped_undefined_sums;
                continue;
            }
            check(true, g, _mm_set1_epi16(extension), h,
                  _mm_set1_epi16((int)sum));
        }
    }
    __m128i r = check(true, _mm_set1_epi16(0), _mm_set1_epi16(1),
                     _mm_set1_epi16(32767), _mm_set1_epi16(32768));
    int values[16]; simd_store(r.words, values);
    assert(values[0] == 65535); // Existing subtraction has no upper clamp.
}

static void test_random_and_chains() {
    for (int n = 0; n < 20000; ++n) {
        int g[16], e[16], h[16], o[16];
        for (int i = 0; i < 16; ++i) {
            g[i] = (int)(random_word() % 196609) - 65536;
            h[i] = (int)(random_word() % 196609) - 65536;
            e[i] = narrow_signed((int)(random_word() & 65535), 16);
            o[i] = narrow_signed((int)(random_word() & 65535), 16);
        }
        check(true, load(g), load(e), load(h), load(o));
        check(false, load(g), load(e), load(h), load(o));
#if SWBWA_KSW_I16_MODE == SWBWA_KSW_I16_INT32_8
        /* Full word wrap is defined by SIMD, but not by scalar_8's C a-b. */
        for (int i = 0; i < 16; ++i) {
            g[i] = signed_word(random_word()); h[i] = signed_word(random_word());
            e[i] = signed_word(random_word()); o[i] = signed_word(random_word());
        }
        check(true, load(g), load(e), load(h), load(o));
        check(false, load(g), load(e), load(h), load(o));
#endif
    }
    for (int penalty : {-32768, -1, 0, 1, 32767}) {
        __m128i gap = _mm_set1_epi16(0);
        __m128i ext = _mm_set1_epi16(penalty);
        for (int n = 0; n < 1000; ++n)
            gap = check(true, gap, ext, _mm_set1_epi16(32767), ext);
    }
}

/* Independent algebra checks for the production XOR select/clamp operations.
 * No substitute implementation is injected into the included header.
 */
static intv16 xor_select(intv16 a, intv16 b, intv16 choose_b) {
    intv16 zero = 0;
    intv16 mask = simd_vsubw(zero, choose_b);
    return simd_vxorw(a, simd_vandw(simd_vxorw(a, b), mask));
}
static void test_instruction_algebra() {
    const int edges[] = {INT_MIN, INT_MIN + 1, -32768, -255, -1, 0, 1, 127,
                         128, 254, 255, 256, 32767, 65535, INT_MAX - 1, INT_MAX};
    for (int n = 0; n < 20000; ++n) {
        int av[16], bv[16], choose[16];
        for (int i = 0; i < 16; ++i) {
            av[i] = n < 16 ? edges[i] : signed_word(random_word());
            bv[i] = n < 16 ? edges[(i + n) % 16] : signed_word(random_word());
            if ((i + n) % 4 == 0) bv[i] = av[i];
            choose[i] = random_word() & 1;
        }
        __m128i a = load(av), b = load(bv), result, reference;
        intv16 zero = 0, upper = 255, choice = load(choose).words;
        result.words = xor_select(a.words, b.words, choice);
        reference.words = swbwa_i16_select(a.words, b.words, choice);
        equal(result, reference, "test-only xor select");
        result.words = xor_select(a.words, b.words, simd_vcmpltw(a.words, b.words));
        equal(result, _mm_max_epu8(a, b), "test-only xor max");
        result.words = xor_select(a.words, b.words, simd_vcmpltw(b.words, a.words));
        equal(result, _mm_min_epu8(a, b), "test-only xor min");
        intv16 sum = simd_vaddw(a.words, b.words);
        result.words = xor_select(upper, sum, simd_vcmpltw(sum, upper));
        equal(result, _mm_adds_epu8(a, b), "test-only xor saturated add");
        intv16 difference = simd_vsubw(a.words, b.words);
        result.words = simd_vandw(difference,
                                 simd_vsubw(zero, simd_vcmpltw(zero, difference)));
        equal(result, _mm_subs_epu8(a, b), "test-only positive-mask subtract");
        proposal_cases += 5;
    }
}

int main() {
    test_shim();
    test_u8_differences();
    test_penalty_narrowing();
    test_random_and_chains();
    test_instruction_algebra();
    std::printf("i16_mode=%d: %llu actual-helper cases, %llu lane checks; "
                "%llu pre-conversion int-overflow sums excluded: PASS\n",
                SWBWA_KSW_I16_MODE, (unsigned long long)cases,
                (unsigned long long)lane_checks,
                (unsigned long long)skipped_undefined_sums);
    std::printf("test-only select/clamp proposals: %llu vector comparisons: PASS\n",
                (unsigned long long)proposal_cases);
}
'''


def main():
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    with tempfile.TemporaryDirectory(prefix="swbwa-ksw-gap-") as temporary:
        path = Path(temporary)
        (path / "simd.h").write_text(SIMD, encoding="ascii")
        (path / "test.cpp").write_text(HARNESS, encoding="ascii")
        for mode, optimization in ((2, "-O2"), (2, "-O3"), (1, "-O2")):
            executable = path / f"gap-{mode}-{optimization[1:]}"
            command = compiler + [
                "-std=c++11", optimization, "-g", "-Wall", "-Wextra", "-Werror",
                "-Wno-psabi", "-fsanitize=undefined", "-fno-sanitize-recover=all",
                "-DSWBWA_KSW_FUSED_GAP_UPDATE=1",
                f"-DSWBWA_KSW_I16_MODE={mode}",
                "-I", str(path), "-I", str(ROOT / "include"),
                "-I", str(ROOT), str(path / "test.cpp"),
                "-o", str(executable),
            ]
            print(f"Compiling actual scalar_sse.h: i16_mode={mode} {optimization}",
                  flush=True)
            subprocess.run(command, check=True)
            subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    main()
