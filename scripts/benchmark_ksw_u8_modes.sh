#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}/.."

RESULT_DIR=${RESULT_DIR:-correctness_results/ksw_u8_modes_bigdata_20260827}
DATA_DIR=${DATA_DIR:-../data/bwa_test_big_data}
REF=${REF:-../data/GRCh38.d1.vd1.fa}
BUILD_JOBS=${BUILD_JOBS:-8}

MODES=(int32_16 float16_16 float16_32)
DATASETS=(ERR1203383 small_SRR7963242 SRR2496709)
SUMMARY="${RESULT_DIR}/summary.tsv"

mkdir -p "$RESULT_DIR"
printf 'mode\tdataset\tstage2_s\tpart3_s\tlines\tbytes\tmd5\treference_md5\tstatus\n' \
    > "$SUMMARY"

expected_md5()
{
    local dataset=$1

    awk -F '\t' -v dataset="$dataset" \
        '$1 == dataset && $2 == "PE" { print $3; exit }' \
        scripts/bigdata_expected_md5.tsv
}

for mode in "${MODES[@]}"; do
    mode_dir="${RESULT_DIR}/${mode}"
    mkdir -p "$mode_dir"

    echo "[BUILD] KSW_U8_MODE=$mode"
    ./build_cross.sh "$BUILD_JOBS" \
        EXEC_MODE=cgs_cross \
        CPE_ALLOCATOR=pool \
        HOST_MALLOC_WRAPPER=1 \
        HOST_MALLOC_STATS=0 \
        CPE_PROFILE=0 \
        KSW_U8_MODE="$mode" \
        USE_MPI=0 \
        > "${mode_dir}/build.log" 2>&1

    for dataset in "${DATASETS[@]}"; do
        r1="${DATA_DIR}/${dataset}_1.fastq"
        r2="${DATA_DIR}/${dataset}_2.fastq"
        sam="${mode_dir}/${dataset}_PE.sam"
        log="${mode_dir}/${dataset}_PE.run.log"
        expected=$(expected_md5 "$dataset")

        [[ -n "$expected" ]] || {
            echo "missing reference MD5 for $dataset PE" >&2
            exit 1
        }
        rm -f "$sam" "$log"

        echo "[RUN] mode=$mode dataset=$dataset PE"
        bsub -I -b -q q_share \
            -n 1 -cgsp 64 \
            -share_size 2000 -mpecg 6 \
            -xmalloc -cross_size 42000 \
            -cache_size 128 -priv_size 16 \
            -J "ksw_${mode}_${dataset}" -o "$log" \
            ./SWBWA mem -v 4 -t 1 -1 \
            -o "$sam" -I 170,80,500,1 \
            "$REF" "$r1" "$r2"

        [[ -f "$sam" ]] || {
            echo "SAM output missing: $sam" >&2
            exit 1
        }

        stage2=$(awk '/stage 2 - align reads/{value=$(NF-1)} END {print value}' "$log")
        part3=$(awk '/part 3 - CPE alignment/{value=$(NF-1)} END {print value}' "$log")
        lines=$(wc -l < "$sam")
        bytes=$(wc -c < "$sam")
        actual=$(md5sum "$sam" | awk '{print $1}')
        status=MISMATCH
        [[ "$actual" == "$expected" ]] && status=PASS

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$mode" "$dataset" "$stage2" "$part3" "$lines" "$bytes" \
            "$actual" "$expected" "$status" | tee -a "$SUMMARY"
        rm -f "$sam"
    done
done

echo "[DONE] $SUMMARY"
