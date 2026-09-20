#!/usr/bin/env bash
# Verify that the profiled (CPE_PROFILE=1) builds reproduce the reference SAM MD5.
#
# Two bsub jobs per variant:
#   1. SWBWA itself, as the DIRECT command (mandatory on Sunway: wrapping it in a
#      shell makes the CPE launcher set up the wrong mapping and it segfaults).
#   2. md5sum on the produced SAM, also as its own bsub job, so the several-GB
#      read stays off the login node.
# Expected MD5 is the header-free SAM reference from
# scripts/bigdata_expected_md5.tsv.
#
# Usage: ./md5_check.sh [cross|single] [dataset] [PE|SE] [expected_md5]
set -euo pipefail

AB_ROOT=${AB_ROOT:?Set AB_ROOT to the experiment directory}
QUEUE=${QUEUE:-q_share}
CFG=${1:-cross}
DS=${2:-ERR1203383}
MODE=${3:-PE}
EXPECT=${4:-}
DATA_ROOT=${DATA_ROOT:?Set DATA_ROOT to the reference and input directory}
[[ "$CFG" == cross || "$CFG" == single ]] || exit 2
[[ "$MODE" == PE || "$MODE" == SE ]] || exit 2
[[ "$EXPECT" =~ ^[0-9a-f]{32}$ ]] || { echo "Expected MD5 is required" >&2; exit 2; }
mkdir -p "$AB_ROOT/logs"
exec 9>"$AB_ROOT/driver.lock"
flock -n 9 || { echo "An A/B driver is already active" >&2; exit 1; }

CROSS_ARGS=(-n 1 -cgsp 64 -share_size 2000 -mpecg 6 -xmalloc -cross_size 42000
            -cache_size 128 -priv_size 16)
SINGLE_ARGS=(-n 1 -cgsp 64 -share_size 12000 -cache_size 128 -priv_size 16)
ALLOC=("${SINGLE_ARGS[@]}")
[[ "$CFG" == cross ]] && ALLOC=("${CROSS_ARGS[@]}")

REF=$DATA_ROOT/GRCh38.d1.vd1.fa
R1=$DATA_ROOT/bwa_test_big_data/${DS}_1.fastq
R2=$DATA_ROOT/bwa_test_big_data/${DS}_2.fastq
mkdir -p "$AB_ROOT/scratch"

echo "### MD5 check: cfg=$CFG ds=$DS mode=$MODE expect=${EXPECT:-<none>}"
overall=0
for variant in old new; do
    rdir="$AB_ROOT/runs/${variant}_${CFG}_p1"
    sam="$AB_ROOT/scratch/md5check_${variant}.sam"
    rlog="$AB_ROOT/logs/md5check_${variant}_run.log"
    mlog="$AB_ROOT/logs/md5check_${variant}_md5.log"
    [[ -x "$rdir/SWBWA" ]] || { echo "[MD5] missing $rdir/SWBWA"; overall=1; continue; }
    [[ ! -e "$sam" ]] || { echo "[MD5] preserving existing $sam" >&2; exit 1; }

    cargs=(mem -v 4 -t 1 -1 -o "$sam" -I 170,80,500,1 "$REF" "$R1")
    [[ "$MODE" == PE ]] && cargs+=("$R2")

    echo "[MD5] run $variant $(date '+%F %T')"
    ( cd "$rdir" && bsub -I -b -q "$QUEUE" "${ALLOC[@]}" \
        -J "md5run_${variant}_${CFG}_${DS}_${MODE}" -o "$rlog" \
        ./SWBWA "${cargs[@]}" ) > "$rlog.out" 2>&1
    echo "[MD5] run $variant rc=$? sam_bytes=$(stat -c%s "$sam" 2>/dev/null || echo 0)"
    grep -q 'has been finished' "$rlog.out" || exit 1
    grep -q '\[main\] CMD:' "$rlog" "$rlog.out" || exit 1

    echo "[MD5] md5sum $variant $(date '+%F %T')"
    bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 -share_size 2000 \
        -cache_size 128 -priv_size 16 \
        -J "md5sum_${variant}_${CFG}_${DS}_${MODE}" -o "$mlog" \
        md5sum "$sam" > "$mlog.out" 2>&1
    grep -q 'has been finished' "$mlog.out" || exit 1

    actual=$(awk '{print $1}' "$mlog.out" | grep -E '^[0-9a-f]{32}$' | head -1)
    if [[ -n "$EXPECT" && "$actual" == "$EXPECT" ]]; then
        echo "[MD5] $variant PASS $actual"
        rm -f "$sam"
        echo "[MD5] $variant sam removed"
    else
        echo "[MD5] $variant MISMATCH expect=$EXPECT actual=${actual:-<none>}"
        overall=1
    fi
done
echo "### MD5 check done overall_rc=$overall"
exit "$overall"
