#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
"${CC:-cc}" -std=c11 -Wall -Wextra -Werror -Wno-unused-parameter \
    -Wno-deprecated-declarations -fsanitize=address,undefined -g -O1 \
    -Itests/stubs -Isrc/slave -Iinclude -include include/swbwa_config.h \
    -DSWBWA_EXEC_MODE=SWBWA_EXEC_SINGLE_CG \
    -DSWBWA_CPE_ALLOC_MODE=SWBWA_CPE_ALLOC_SYSTEM \
    tests/test_cpe_allocator.c src/slave/malloc_wrap.c -lm -o "$tmp/allocator"
"$tmp/allocator"
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_discard_hash.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_review_tools.py
