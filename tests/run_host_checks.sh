#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
"${CC:-cc}" -std=c11 -Wall -Wextra -Werror -Wno-unused-parameter \
    -Wno-deprecated-declarations -fsanitize=address,undefined -fno-sanitize-recover=all -g -O1 \
    -Itests/stubs -Isrc/slave -Iinclude -include include/swbwa_config.h \
    -DSWBWA_EXEC_MODE=SWBWA_EXEC_SINGLE_CG \
    -DSWBWA_CPE_ALLOC_MODE=SWBWA_CPE_ALLOC_SYSTEM \
    tests/test_cpe_allocator.c src/slave/malloc_wrap.c -lm -o "$tmp/allocator"
"$tmp/allocator"
"${CC:-cc}" -std=gnu11 -Wall -Wextra -Werror -O2 -pthread \
    -Iinclude -DSWBWA_HOST_MPE_THREADS=6 tests/test_host_workers.c \
    src/host/swbwa_host_workers.c -o "$tmp/host_workers"
"$tmp/host_workers"
for mpi in 0 1; do
    "${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -Wno-unused-parameter -pthread -O1 -g \
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-pie -no-pie \
        -Iinclude -include include/swbwa_config.h -DSWBWA_USE_MPI="$mpi" \
        -DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_DYNAMIC \
        tests/test_pipeline_queue.cpp src/host/kthread.cpp -o "$tmp/pipeline_queue"
    "$tmp/pipeline_queue"
done
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_discard_hash.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_output_discard.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_cpe_discard_digest.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_optimization_defaults.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_review_tools.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_cross_stack.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_parallel_input.py
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_ordered_sam_check.py
