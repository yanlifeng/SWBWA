#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONDONTWRITEBYTECODE=1
for test in \
    test_ldm_modes.py \
    test_ldm_allocator.py \
    test_ksw_fused_gap.py \
    test_ksw_fused_gap_native.py \
    test_ksw_xor_select_native.py \
    test_ksw_u8_native.py \
    test_ksw_i16_native.py \
    test_ksw_extend2_qp_native.py \
    test_ksw_extend2_ldm_native.py \
    test_chain_reuse_smem_ldm.py \
    test_ksw_modes_compile.py
do
    printf '\n[CHECK] %s\n' "$test"
    python3 "tests/$test"
done
