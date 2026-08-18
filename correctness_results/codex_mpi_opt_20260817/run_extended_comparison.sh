#!/usr/bin/env bash
set -euo pipefail

if (($# != 3)); then
    echo "usage: $0 VERSION SOURCE_DIR OUTPUT_MODE" >&2
    exit 2
fi

VERSION=$1
SOURCE_DIR=$2
OUTPUT_MODE=$3

case "$OUTPUT_MODE" in
    split|single_unordered) ;;
    *) echo "invalid output mode: $OUTPUT_MODE" >&2; exit 2 ;;
esac

if [[ "$SOURCE_DIR" == "." ]]; then
    DATA=../data
    REF=../data/GRCh38.d1.vd1.fa
    SORTER=../fucking_sam_sort
    RESULT_BASE=correctness_results/codex_mpi_opt_20260817/extended_comparison
else
    DATA=../../data
    REF=../../data/GRCh38.d1.vd1.fa
    SORTER=../../fucking_sam_sort
    RESULT_BASE=../correctness_results/codex_mpi_opt_20260817/extended_comparison
fi

OUTPUT_DIR=${RESULT_BASE}/${VERSION}/${OUTPUT_MODE}
SUBSET_DIR=${OUTPUT_DIR}/subsets
FULL_DIR=${OUTPUT_DIR}/full

cd "$SOURCE_DIR"
mkdir -p "$SUBSET_DIR" "$FULL_DIR"

echo "[BATCH] version=$VERSION output=$OUTPUT_MODE source=$SOURCE_DIR"

DATA=$DATA REF=$REF SORTER=$SORTER VERBOSE=4 \
    bash scripts/correctness.sh run \
        --profile mpi \
        --io-mode no1 \
        --datasets SRR2496709,ERR1203383 \
        --mpi-input dynamic \
        --output-mode "$OUTPUT_MODE" \
        --output-dir "$SUBSET_DIR"

DATA=$DATA REF=$REF SORTER=$SORTER VERBOSE=4 \
    bash scripts/correctness.sh verify \
        --datasets SRR2496709,ERR1203383 \
        --mpi-input dynamic \
        --output-mode "$OUTPUT_MODE" \
        --output-dir "$SUBSET_DIR"

FULL_BASE=${FULL_DIR}/SRR2496709_FULL_PE.sam
FULL_COMBINED=${FULL_DIR}/combined_SRR2496709_FULL_PE.sam
FULL_NORMALIZED=${FULL_DIR}/normalized_SRR2496709_FULL_PE.sam
FULL_RUN_LOG=${FULL_DIR}/SRR2496709_FULL_PE.run.log
FULL_SORT_LOG=${FULL_DIR}/SRR2496709_FULL_PE.sort.log
FULL_MD5=${FULL_DIR}/md5_results.tsv

rm -f "$FULL_BASE" "$FULL_COMBINED" "$FULL_NORMALIZED" \
      "$FULL_RUN_LOG" "$FULL_SORT_LOG" "$FULL_MD5" \
      "${FULL_BASE%.sam}".rank*.sam

bsub -I -b -q q_share -N 1 -np 6 -cgsp 64 \
    -J "full_${VERSION}_${OUTPUT_MODE}" \
    -share_size 12000 -cache_size 128 -priv_size 16 \
    -o "$FULL_RUN_LOG" \
    ./SWBWA mem -v 4 -t 1 \
        -o "$FULL_BASE" -I 170,80,500,1 \
        "$REF" "$DATA/SRR2496709_1.fastq" "$DATA/SRR2496709_2.fastq"

if [[ "$OUTPUT_MODE" == split ]]; then
    rank_files=()
    for rank in 0 1 2 3 4 5; do
        rank_file=$(printf '%s.rank%06d.sam' "${FULL_BASE%.sam}" "$rank")
        [[ -f "$rank_file" ]] || {
            echo "missing split output: $rank_file" >&2
            exit 1
        }
        rank_files+=("$rank_file")
    done
    cat "${rank_files[@]}" > "$FULL_COMBINED"
    rm -f "${rank_files[@]}"
    FULL_SOURCE=$FULL_COMBINED
else
    [[ -f "$FULL_BASE" ]] || {
        echo "missing single output: $FULL_BASE" >&2
        exit 1
    }
    FULL_SOURCE=$FULL_BASE
fi

bsub -I -b -q q_share -n 1 -cgsp 64 \
    -J "sort_full_${VERSION}_${OUTPUT_MODE}" \
    -share_size 14000 -cache_size 128 -priv_size 16 \
    -o "$FULL_SORT_LOG" \
    "$SORTER" "$FULL_SOURCE" "$FULL_NORMALIZED"

[[ -s "$FULL_NORMALIZED" ]] || {
    echo "normalized full SAM was not generated" >&2
    exit 1
}

MD5=$(md5sum "$FULL_NORMALIZED" | awk '{print $1}')
printf 'version\toutput_mode\tdataset\tread_mode\tmd5\n' > "$FULL_MD5"
printf '%s\t%s\tSRR2496709_FULL\tPE\t%s\n' \
       "$VERSION" "$OUTPUT_MODE" "$MD5" >> "$FULL_MD5"
echo "[FULL MD5] version=$VERSION output=$OUTPUT_MODE md5=$MD5"

rm -f "$FULL_BASE" "$FULL_COMBINED" "$FULL_NORMALIZED" \
      "${FULL_BASE%.sam}".rank*.sam

echo "[BATCH PASS] version=$VERSION output=$OUTPUT_MODE"
