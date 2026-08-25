#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}/../.."

RESULT_DIR=correctness_results/chunk_size_single_system_20260824
REF=../data/GRCh38.d1.vd1.fa
DATA=../data
SIZES_MIB=(1 2 4 8 16 32 64)
DATASETS=(SRR7963242 SRR2496709 ERR1203383)

mkdir -p "$RESULT_DIR"

make clean
make -j8 USE_MPI=0 EXEC_MODE=single CPE_ALLOCATOR=system \
    HOST_MALLOC_WRAPPER=1 HOST_MALLOC_STATS=0 CPE_PROFILE=0

for dataset in "${DATASETS[@]}"; do
    input="${DATA}/${dataset}_1.fastq_1"
    [[ -f "$input" ]] || {
        echo "missing input: $input" >&2
        exit 1
    }

    for size_mib in "${SIZES_MIB[@]}"; do
        chunk_bytes=$((size_mib * 1024 * 1024))
        case_dir="${RESULT_DIR}/${dataset}"
        log="${case_dir}/chunk_${size_mib}mib.log"
        marker="${case_dir}/chunk_${size_mib}mib.ok"

        mkdir -p "$case_dir"
        if [[ -f "$marker" ]]; then
            echo "[SKIP] ${dataset} K=${size_mib} MiB"
            continue
        fi

        echo "[RUN ] ${dataset} K=${size_mib} MiB"
        bsub -I -b -q q_share -n 1 -cgsp 64 \
            -J "k${size_mib}_${dataset}" \
            -share_size 12000 -cache_size 128 -priv_size 16 \
            -o "$log" \
            ./SWBWA mem -v 4 -t 1 -1 -K "$chunk_bytes" \
            -o /dev/null -I 170,80,500,1 "$REF" "$input"

        grep -q "SWBWA Timing Report" "$log"
        touch "$marker"
    done
done

touch "${RESULT_DIR}/run.ok"
