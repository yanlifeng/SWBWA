"""Compile/type-check actual KSW alignment and extension in all three u8 modes.

Run: python3 tests/test_ksw_modes_compile.py. Checks both i16 modes, dual-forward
off/on, and all fused/XOR combinations. The full pre-global-alignment ksw.c
section is compiled unchanged, including the real header, forward/reverse
wrappers, and scalar extension. Verify default-off LDM and enabled-LDM builds.
Also check umbrella defaults and explicit per-feature overrides in every mode.

Local GCC 11 lacks _Float16. A declaration-only, two-byte half adapter permits
object compilation; these objects must not be linked or executed. This checks
conditional/type compatibility, NOT Sunway intrinsic lowering or FP behavior.
No performance experiment is performed. All artifacts use a temporary directory.
"""

import hashlib
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from test_ksw_i16_native import ROOT, SIMD


HALF_DECLARATIONS = r'''
#ifndef SWBWA_TEST_HALF_DECLARATIONS_H
#define SWBWA_TEST_HALF_DECLARATIONS_H
// Declaration-only adapter: preserve storage width and numeric conversions.
// Intentionally no definitions; an accidentally linked FP test must fail.
struct swbwa_test_half {
    uint16_t bits;
    swbwa_test_half() = default;
    swbwa_test_half(float);
    operator float() const;
};
typedef swbwa_test_half _Float16;
struct alignas(64) float16v32 {
    _Float16 lanes[32];
    float16v32() = default;
    float16v32(_Float16);
    float16v32(float);
    float16v32(int);
};
static_assert(sizeof(_Float16) == 2, "half adapter must be two bytes");
static_assert(sizeof(float16v32) == 64, "half vector must be 64 bytes");
void simd_load(float16v32 &, const _Float16 *);
void simd_store(float16v32, _Float16 *);
float16v32 simd_vaddh(float16v32, float16v32);
float16v32 simd_vsubh(float16v32, float16v32);
float16v32 simd_smaxh(float16v32, float16v32);
float16v32 simd_sminh(float16v32, float16v32);
_Float16 simd_reduc_smaxh(float16v32);
#endif
'''

PREAMBLE = r'''
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "swbwa_config.h"
#include "simd.h"
#undef __SSE2__
#undef __ARM_NEON
'''

ASSERTIONS = r'''
static_assert(sizeof(__m128i) == 64, "all modes retain full physical vectors");
static_assert(SWBWA_KSW_U8_LANES == SWBWA_KSW_U8_LOGICAL_LANES,
              "query initialization and kernel stripes must agree");
static_assert(SWBWA_KSW_U8_LANES == (SWBWA_KSW_U8_MODE == 3 ? 32 : 16),
              "the toggles must not change logical lane count");
static_assert(SWBWA_KSW_EXTEND2_LDM_MAX_BYTES == SWBWA_TEST_EXPECT_LDM,
              "extension LDM must default off and honor explicit overrides");
#ifdef SWBWA_TEST_EXPECT_FEATURES
static_assert(SWBWA_KSW_FUSED_GAP_UPDATE == SWBWA_TEST_EXPECT_FEATURES,
              "fused gap must honor umbrella defaults and explicit overrides");
static_assert(SWBWA_KSW_XOR_SELECT == SWBWA_TEST_EXPECT_FEATURES,
              "XOR must honor umbrella defaults and explicit overrides");
static_assert(SWBWA_KSW_EXTEND2_QP_UNROLL == SWBWA_TEST_EXPECT_FEATURES,
              "QP unroll must honor umbrella defaults and explicit overrides");
#endif
'''

FP_U8_HELPERS = (
    "_mm_set1_epi32", "_mm_load_si128", "_mm_store_si128", "m128i_allzero",
    "_mm_slli_si128", "_mm_max_epu8", "_mm_min_epu8", "m128i_max_u8",
    "_mm_set1_epi8", "_mm_adds_epu8", "_mm_subs_epu8",
)


def body(preprocessed, name):
    """Find a function definition, not a prototype, in comment-free CPP output."""
    matches = []
    for match in re.finditer(r"\b" + re.escape(name) + r"\s*\(", preprocessed):
        depth = 1
        after = match.end()
        while depth:
            if preprocessed[after] == "(":
                depth += 1
            elif preprocessed[after] == ")":
                depth -= 1
            after += 1
        while after < len(preprocessed) and preprocessed[after].isspace():
            after += 1
        if after < len(preprocessed) and preprocessed[after] == "{":
            matches.append((match.start(), after))
    assert len(matches) == 1, (name, len(matches))
    begin, opening = matches[0]
    depth = 1
    end = opening + 1
    while depth:
        if preprocessed[end] == "{":
            depth += 1
        elif preprocessed[end] == "}":
            depth -= 1
        end += 1
    return preprocessed[begin:end]


def main():
    assert body("int f(int x); int f(int x) { return x; } "
                "if (f(g(x))) { use(x); }", "f") == "f(int x) { return x; }"
    source = (ROOT / "src/slave/ksw.c").read_text()
    marker = "\n/********************\n * Global alignment *"
    assert source.count(marker) == 1
    # Exclude only the unrelated global-alignment section with C-only casts.
    source = source[:source.index(marker)]
    compiler = shlex.split(os.environ.get("CXX", "g++"))
    count = 0
    with tempfile.TemporaryDirectory(prefix="swbwa-ksw-modes-") as temporary:
        path = Path(temporary)
        (path / "simd.h").write_text(SIMD + HALF_DECLARATIONS, encoding="ascii")
        (path / "test.cpp").write_text(
            PREAMBLE + '\n#line 1 "src/slave/ksw.c"\n' + source + ASSERTIONS,
            encoding="ascii")
        base = compiler + [
            "-std=c++11", "-O0", "-Wall", "-Wextra", "-Werror", "-Wno-psabi",
            "-DSWBWA_ENABLE_CPE_PROFILE=0", "-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=0",
            "-I", str(path), "-I", str(ROOT / "src/slave"),
            "-I", str(ROOT / "include"), "-I", str(ROOT / "tests/stubs"),
        ]
        print("Compile-only: declaration-only half adapter; no FP execution",
              flush=True)
        for mode in (1, 2, 3):
            for i16 in (1, 2):
                for dual in (0, 1):
                    baseline = None
                    for fused, xor in ((0, 0), (1, 0), (0, 1), (1, 1)):
                        name = f"u8-{mode}-i16-{i16}-dual-{dual}-f{fused}-x{xor}"
                        flags = [f"-DSWBWA_KSW_U8_MODE={mode}",
                                 f"-DSWBWA_KSW_I16_MODE={i16}",
                                 f"-DSWBWA_ENABLE_MATESW_DUAL_FORWARD={dual}",
                                 f"-DSWBWA_KSW_FUSED_GAP_UPDATE={fused}",
                                 f"-DSWBWA_KSW_XOR_SELECT={xor}",
                                 "-DSWBWA_TEST_EXPECT_LDM=0"]
                        subprocess.run(base + flags + [
                            "-c", str(path / "test.cpp"), "-o", str(path / "test.o")
                        ], check=True)
                        pp = subprocess.check_output(base + flags + [
                            "-E", "-P", str(path / "test.cpp")
                        ], text=True)
                        functions = {n: body(pp, n) for n in
                                     ("ksw_u8", "ksw_i16", "ksw_qinit_impl",
                                      "ksw_align2_matesw_dual_forward", "ksw_extend2")}
                        functions["dual"] = (body(pp, "ksw_u8_dual_forward")
                                             if mode == 1 and dual else "")
                        assert ("swbwa_ksw_gap_update_words(" in functions["ksw_u8"]
                                ) == bool(mode == 1 and fused), name
                        assert ("swbwa_ksw_i16_gap_update(" in functions["ksw_i16"]
                                ) == bool(i16 == 2 and fused), name
                        assert "swbwa_ldm_alloc(" not in functions["ksw_extend2"], name
                        assert ("swbwa_ksw_xor_select_words(" in body(pp, "swbwa_i16_select")
                                ) == bool(xor), name
                        if fused:
                            for n in ("swbwa_ksw_gap_update_words",
                                      "swbwa_ksw_i16_gap_update"):
                                definition = body(pp, n)
                                assert ".words" in definition and ".val" not in definition, name
                        if mode != 1:
                            functions.update({n: body(pp, n) for n in FP_U8_HELPERS})
                        if baseline is None:
                            baseline = functions
                        else:
                            for n in ("ksw_qinit_impl", "dual",
                                      "ksw_align2_matesw_dual_forward", "ksw_extend2"):
                                assert functions[n] == baseline[n], (name, n)
                            if mode != 1:
                                for n in ("ksw_u8",) + FP_U8_HELPERS:
                                    assert functions[n] == baseline[n], (name, n)
                        count += 1
                    print(f"u8={mode} i16={i16} dual={dual}: four toggles compile; "
                          "conditional checks PASS", flush=True)
        for mode in (1, 2, 3):
            for qp in (0, 1):
                flags = [f"-DSWBWA_KSW_U8_MODE={mode}", "-DSWBWA_KSW_I16_MODE=2",
                         "-DSWBWA_ENABLE_MATESW_DUAL_FORWARD=1",
                         "-DSWBWA_KSW_FUSED_GAP_UPDATE=1", "-DSWBWA_KSW_XOR_SELECT=1",
                         "-DSWBWA_KSW_EXTEND2_LDM_MAX_BYTES=8192",
                         "-DSWBWA_TEST_EXPECT_LDM=8192",
                         f"-DSWBWA_KSW_EXTEND2_QP_UNROLL={qp}"]
                subprocess.run(base + flags + [
                    "-c", str(path / "test.cpp"), "-o", str(path / "test.o")
                ], check=True)
                pp = subprocess.check_output(base + flags + [
                    "-E", "-P", str(path / "test.cpp")
                ], text=True)
                extension = body(pp, "ksw_extend2")
                assert "swbwa_ldm_alloc(" in extension
                assert "swbwa_ldm_release(" in extension
                count += 1
            print(f"u8={mode}: enabled LDM=8192, QP off/on compile PASS", flush=True)
        for mode in (1, 2, 3):
            for label, umbrella, override in (
                ("standalone", None, None), ("umbrella-off", 0, None),
                ("umbrella-on", 1, None), ("override-off", 1, 0),
                ("override-on", 0, 1),
            ):
                inherited = umbrella or 0
                expected = inherited if override is None else override
                cap = 4096 * inherited if override is None else 8192 * override
                flags = [f"-DSWBWA_KSW_U8_MODE={mode}", "-DSWBWA_KSW_I16_MODE=2",
                         "-DSWBWA_ENABLE_MATESW_DUAL_FORWARD=1",
                         f"-DSWBWA_TEST_EXPECT_LDM={cap}",
                         f"-DSWBWA_TEST_EXPECT_FEATURES={expected}"]
                if umbrella is not None:
                    flags.append(f"-DSWBWA_ENABLE_CPE_KERNEL_OPT={umbrella}")
                if override is not None:
                    flags += [f"-DSWBWA_KSW_FUSED_GAP_UPDATE={override}",
                              f"-DSWBWA_KSW_XOR_SELECT={override}",
                              f"-DSWBWA_KSW_EXTEND2_QP_UNROLL={override}",
                              f"-DSWBWA_KSW_EXTEND2_LDM_MAX_BYTES={cap}"]
                subprocess.run(base + flags + [
                    "-c", str(path / "test.cpp"), "-o", str(path / "test.o")
                ], check=True)
                pp = subprocess.check_output(base + flags + [
                    "-E", "-P", str(path / "test.cpp")
                ], text=True)
                assert ("swbwa_ksw_gap_update_words(" in body(pp, "ksw_u8")
                        ) == bool(mode == 1 and expected), (mode, label)
                assert ("swbwa_ksw_i16_gap_update(" in body(pp, "ksw_i16")
                        ) == bool(expected), (mode, label)
                assert ("swbwa_ksw_xor_select_words(" in body(pp, "swbwa_i16_select")
                        ) == bool(expected), (mode, label)
                assert ("swbwa_ldm_alloc(" in body(pp, "ksw_extend2")
                        ) == bool(cap), (mode, label)
                count += 1
            print(f"u8={mode}: standalone/umbrella defaults and overrides PASS", flush=True)
        print(f"PASS: {count} object compiles; alignment/extension-section SHA256="
              f"{hashlib.sha256(source.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
