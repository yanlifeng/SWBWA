"""Compare actual scalar_sse.h with XOR disabled/enabled in one host process.

Uses the GCC-vector SIMD shim from test_ksw_fused_gap_native.py. No source
rewriting or substitute production helpers; both header versions are checked
against independent scalar oracles. Sunway lowering/performance is not tested.
"""

import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from test_ksw_fused_gap_native import SIMD


ROOT = Path(__file__).resolve().parents[1]
HARNESS = r'''
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include "swbwa_config.h"
#include "simd.h"

namespace baseline {
#include "src/slave/scalar_sse.h"
static_assert(SWBWA_KSW_XOR_SELECT == 0, "XOR must default off");
}
#undef SCALAR_SSE_H
#undef SWBWA_KSW_XOR_SELECT
#define SWBWA_KSW_XOR_SELECT 1
namespace candidate {
#include "src/slave/scalar_sse.h"
}

static uint64_t comparisons, lanes;
static uint32_t state = 7963242;

static uint32_t random_word() {
    state ^= state << 13; state ^= state >> 17; state ^= state << 5;
    return state;
}
static int32_t signed_word(uint32_t bits) {
    int32_t r;
    std::memcpy(&r, &bits, sizeof(r));
    return r;
}
static int32_t add_word(int32_t a, int32_t b) {
    return signed_word((uint32_t)a + (uint32_t)b);
}
static int32_t sub_word(int32_t a, int32_t b) {
    return signed_word((uint32_t)a - (uint32_t)b);
}
static int narrow16(int n) { return (n & 65535) < 32768 ? (n & 65535) : (n & 65535) - 65536; }

enum operation { MAXIMUM, MINIMUM, ADD_U8, SUB_ZERO, SELECT, ADD_I16 };
static intv16 oracle(operation op, intv16 a, intv16 b, intv16 choose, int active = 16) {
    intv16 r = 0;
    for (int i = 0; i < active; ++i) {
        int x = a.lanes[i], y = b.lanes[i];
        switch (op) {
        case MAXIMUM: r.lanes[i] = std::max(x, y); break;
        case MINIMUM: r.lanes[i] = std::min(x, y); break;
        case ADD_U8: r.lanes[i] = std::min(add_word(x, y), 255); break;
        case SUB_ZERO: r.lanes[i] = std::max(sub_word(x, y), 0); break;
        case SELECT: r.lanes[i] = choose.lanes[i] ? y : x; break;
        case ADD_I16: r.lanes[i] = std::max(-32768, std::min(x + y, 32767)); break;
        }
    }
    return r;
}
static void verify(const char *name, intv16 old, intv16 now, intv16 expected) {
    ++comparisons;
    for (int i = 0; i < 16; ++i) {
        ++lanes;
        if (old.lanes[i] != now.lanes[i] || now.lanes[i] != expected.lanes[i]) {
            std::fprintf(stderr, "%s comparison=%llu lane=%d: off=%d on=%d expected=%d\n",
                         name, (unsigned long long)comparisons, i,
                         old.lanes[i], now.lanes[i], expected.lanes[i]);
            std::abort();
        }
    }
}

static void check_words(intv16 av, intv16 bv, intv16 choose) {
    baseline::__m128i a, b;
    candidate::__m128i x, y;
    a.words = x.words = av; b.words = y.words = bv;
    intv16 maximum = oracle(MAXIMUM, av, bv, choose);
    intv16 minimum = oracle(MINIMUM, av, bv, choose);
    verify("max words", baseline::_mm_max_intv16(av, bv),
           candidate::_mm_max_intv16(av, bv), maximum);
    verify("min words", baseline::_mm_min_intv16(av, bv),
           candidate::_mm_min_intv16(av, bv), minimum);
    verify("max u8", baseline::_mm_max_epu8(a, b).words,
           candidate::_mm_max_epu8(x, y).words, maximum);
    verify("min u8", baseline::_mm_min_epu8(a, b).words,
           candidate::_mm_min_epu8(x, y).words, minimum);
    verify("add u8", baseline::_mm_adds_epu8(a, b).words,
           candidate::_mm_adds_epu8(x, y).words, oracle(ADD_U8, av, bv, choose));
    verify("sub u8", baseline::_mm_subs_epu8(a, b).words,
           candidate::_mm_subs_epu8(x, y).words, oracle(SUB_ZERO, av, bv, choose));
    verify("i16 select", baseline::swbwa_i16_select(av, bv, choose),
           candidate::swbwa_i16_select(av, bv, choose), oracle(SELECT, av, bv, choose));
    verify("i16 max words", baseline::swbwa_i16_max_words(av, bv),
           candidate::swbwa_i16_max_words(av, bv), maximum);
    verify("i16 min words", baseline::swbwa_i16_min_words(av, bv),
           candidate::swbwa_i16_min_words(av, bv), minimum);
    verify("i16 masked max", baseline::_mm_max_epi16(a, b).words,
           candidate::_mm_max_epi16(x, y).words, oracle(MAXIMUM, av, bv, choose, 8));
}

static void check_gaps(intv16 gv, intv16 ev, intv16 hv, intv16 ov) {
    baseline::__m128i g, e, h, o;
    candidate::__m128i xg, xe, xh, xo;
    g.words = xg.words = gv; e.words = xe.words = ev;
    h.words = xh.words = hv; o.words = xo.words = ov;
    intv16 expected = 0;
    for (int i = 0; i < 16; ++i)
        expected.lanes[i] = std::max(0, std::max(sub_word(gv.lanes[i], ev.lanes[i]),
                                                sub_word(hv.lanes[i], ov.lanes[i])));
    verify("fused u8", baseline::swbwa_ksw_gap_update_words(g, e, h, o).words,
           candidate::swbwa_ksw_gap_update_words(xg, xe, xh, xo).words, expected);
    for (int i = 8; i < 16; ++i) expected.lanes[i] = 0;
    verify("fused i16", baseline::swbwa_ksw_i16_gap_update(g, e, h, o).words,
           candidate::swbwa_ksw_i16_gap_update(xg, xe, xh, xo).words, expected);
    verify("i16 unsigned subtraction", baseline::_mm_subs_epu16(g, e).words,
           candidate::_mm_subs_epu16(xg, xe).words, oracle(SUB_ZERO, gv, ev, 0, 8));
}

static void test_u8_and_word_edges() {
    intv16 a, b, choose;
    for (int x = 0; x < 256; ++x) {
        for (int y = 0; y < 256; ++y) {
            for (int i = 0; i < 16; ++i) {
                a.lanes[i] = (x + i) & 255;
                b.lanes[i] = (y + 3 * i) & 255;
                choose.lanes[i] = (x + y + i) & 1;
            }
            check_words(a, b, choose);
        }
    }
    const int edges[] = {INT_MIN, INT_MIN + 1, -32768, -255, -1, 0, 1, 127,
                         128, 254, 255, 256, 32767, 65535, INT_MAX - 1, INT_MAX};
    for (int n = 0; n < 20000; ++n) {
        for (int i = 0; i < 16; ++i) {
            a.lanes[i] = n < 16 ? edges[i] : signed_word(random_word());
            b.lanes[i] = n < 16 ? edges[(i + n) % 16] : signed_word(random_word());
            if ((i + n) % 4 == 0) b.lanes[i] = a.lanes[i];
            choose.lanes[i] = random_word() & 1;
        }
        check_words(a, b, choose);
    }
}

static void test_i16_narrowing_and_clamps() {
    const int values[] = {-32768, -1, 0, 1, 255, 32767, 32768, 65535,
                         17, -3, 90000, -65535, 0, 1, 32767, -32768};
    intv16 g, h, a, b;
    for (int i = 0; i < 16; ++i) { g.lanes[i] = values[i]; h.lanes[i] = values[15-i]; }
    for (int raw = 0; raw < 65536; ++raw) {
        baseline::__m128i e = baseline::_mm_set1_epi16(raw);
        candidate::__m128i xe = candidate::_mm_set1_epi16(raw);
        intv16 expected = 0;
        for (int i = 0; i < 8; ++i) expected.lanes[i] = narrow16(raw);
        verify("penalty narrowing", e.words, xe.words, expected);
        check_gaps(g, e.words, h, baseline::_mm_set1_epi16(65535 - raw).words);
        for (int i = 0; i < 16; ++i) {
            a.lanes[i] = narrow16(raw + i * 257);
            b.lanes[i] = narrow16(65535 - raw + i * 129);
        }
        baseline::__m128i ba, bb;
        candidate::__m128i xa, xb;
        ba.words = xa.words = a; bb.words = xb.words = b;
        verify("signed i16 add clamp", baseline::_mm_adds_epi16(ba, bb).words,
               candidate::_mm_adds_epi16(xa, xb).words, oracle(ADD_I16, a, b, 0, 8));
    }
    /* Do not reinterpret subtraction results as uint16_t or clamp at 32767. */
    auto r = candidate::_mm_subs_epu16(candidate::_mm_set1_epi16(32767),
                                      candidate::_mm_set1_epi16(32768));
    assert(r.words.lanes[0] == 65535);
}

int main() {
    test_u8_and_word_edges();
    test_i16_narrowing_and_clamps();
    std::printf("i16_mode=%d: %llu off/on/oracle vector comparisons, %llu lanes: PASS\n",
                SWBWA_KSW_I16_MODE, (unsigned long long)comparisons,
                (unsigned long long)lanes);
}
'''


def main():
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    with tempfile.TemporaryDirectory(prefix="swbwa-xor-select-") as temporary:
        path = Path(temporary)
        (path / "simd.h").write_text(SIMD, encoding="ascii")
        source = path / "test.cpp"
        source.write_text(HARNESS, encoding="ascii")
        for mode, optimization in ((2, "-O2"), (2, "-O3"), (1, "-O2")):
            executable = path / f"xor-{mode}-{optimization[1:]}"
            print(f"Actual header XOR off/on: i16_mode={mode} {optimization}", flush=True)
            subprocess.run(compiler + [
                "-std=c++11", optimization, "-g", "-Wall", "-Wextra", "-Werror",
                "-Wno-psabi", "-fsanitize=address,undefined",
                "-fno-sanitize-recover=all", "-fno-omit-frame-pointer",
                "-fno-pie", "-no-pie", "-DSWBWA_KSW_FUSED_GAP_UPDATE=1",
                f"-DSWBWA_KSW_I16_MODE={mode}", "-I", str(path),
                "-I", str(ROOT), "-I", str(ROOT / "include"),
                str(source), "-o", str(executable),
            ], check=True)
            subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    main()
