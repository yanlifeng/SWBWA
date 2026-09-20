#!/usr/bin/env bash
# T9 performance comparison: 4-node bigshare (4 procs, A) vs 4-node single-CG
# (24 procs, B), discard mode.
#
# Why discard is required: stage3 output (single-file pwrite + RMA offset
# reservation) differs by a 6x process count between the two configs, and the
# noise would drown out the stage2 compute difference. discard writes nothing
# to disk and only emits a 64-bit fingerprint.
#
# Interleaved order A B B A A B cancels the monotonic drift of cluster state.
# Log names must match the LOG_RE in scripts/analyze_mpi_discard_profile.py:
#   (.+)_(PE|SE)\.run(\d+)\.log$
# hence <CFG>_big4_SRR7963242_PE.run<NN>.log, with the dataset group parsed as
# <CFG>_big4_SRR7963242.
#
# Run on the Sunway login node (the script invokes bsub itself); it must live
# under a directory on the online1 filesystem.
set -euo pipefail

W=${WORK_ROOT:?Set WORK_ROOT to the prepared experiment directory}
DATA_ROOT=${DATA_ROOT:?Set DATA_ROOT to the reference and input directory}
NODES=${NODES:-1}
[[ "$NODES" =~ ^[1-9][0-9]*$ ]] || exit 2
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REF=$DATA_ROOT/GRCh38.d1.vd1.fa
D=$DATA_ROOT/bwa_test_big_data
R1=$D/big4_SRR7963242_1.fastq
R2=$D/big4_SRR7963242_2.fastq
QUEUE=${QUEUE:-q_share}
TAG=${TAG:-big4_SRR7963242}

LOGDIR=$W/logs/perf
mkdir -p "$LOGDIR"
exec 9>"$W/driver.lock"
flock -n 9 || { echo "An experiment driver is already active" >&2; exit 1; }
OUT=$LOGDIR/perf_results.tsv
[[ ! -e "$OUT" ]] || { echo "Preserving existing results: $OUT" >&2; exit 1; }
printf 'seq\tconfig\trep\tjobid\twall_s\tlog\tjob_run_time\n' > "$OUT"

for f in "$REF" "$R1" "$R2"; do
    [[ -f "$f" ]] || { echo "缺少输入: $f" >&2; exit 1; }
done
for cfg in A B; do
    [[ -x "$W/runs/${cfg}_discard/SWBWA" ]] || { echo "缺少二进制: runs/${cfg}_discard/SWBWA" >&2; exit 1; }
done

# Interleaved order: A B B A A B.
ORDER=(A B B A A B)
declare -A REP_COUNTER=( [A]=0 [B]=0 )

run_one() {
    local cfg="$1" seq="$2"
    local rdir args jname log hint jid wall submit_epoch done_epoch

    REP_COUNTER[$cfg]=$(( ${REP_COUNTER[$cfg]} + 1 ))
    local rep
    rep=$(printf '%02d' "${REP_COUNTER[$cfg]}")
    jname="disc_${cfg}_${rep}"

    if [ "$cfg" = "A" ]; then
        rdir="$W/runs/A_discard"
        args=(-N "$NODES" -np 1 -cgsp 64 -mpecg 6 -share_size 2000 -xmalloc -cross_size 42000
              -cache_size 128 -priv_size 16)
    else
        rdir="$W/runs/B_discard"
        args=(-N "$NODES" -np 6 -cgsp 64 -share_size 12000 -cache_size 128 -priv_size 16)
    fi

    log="$LOGDIR/${cfg}_${TAG}_PE.run${rep}.log"
    hint="$W/results/discard_${cfg}_${rep}.sam"
    [[ ! -e "$hint" ]] || { echo "Preserving existing $hint" >&2; return 1; }

    echo "[T9] seq=$seq/$(( ${#ORDER[@]} )) cfg=$cfg rep=$rep -> $(basename "$log")  $(date '+%T')"
    submit_epoch=$(date +%s)
    local bsub_out
    bsub_out=$(cd "$rdir" && bsub -I -b -q "$QUEUE" "${args[@]}" \
        -J "$jname" -o "$log" \
        ./SWBWA mem -v 4 -t 1 -1 -o "$hint" -I 170,80,500,1 \
        "$REF" "$R1" "$R2" 2>&1) || {
        printf '%s\n' "$bsub_out" >&2
        echo "Submission failed or detached; inspect bjobs before restarting" >&2
        return 1
    }
    printf '%s\n' "$bsub_out" > "$log.submit.log"
    grep -q 'has been finished' "$log.submit.log" || return 1
    done_epoch=$(date +%s)
    wall=$((done_epoch - submit_epoch))
    jid=$(printf '%s\n' "$bsub_out" | sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' | head -1)

    # discard mode must not produce an output file.
    if [[ -e "$hint" ]]; then
        echo "  !! ERROR: discard 模式产生了输出文件 $hint" >&2
        return 1
    else
        echo "  OK: 未产生输出文件"
    fi

    # fingerprint completeness.
    if python3 "$SCRIPT_DIR/../check_discard_hash.py" "$log"; then
        echo "  OK: 指纹完整"
    else
        echo "  !! ERROR: 指纹未通过完整性检查，停止后续作业" >&2
        return 1
    fi

    local jrt
    jrt=$(bjobs -a -noheader -o "run_time" "$jid" 2>/dev/null | head -1 | tr -d ' ' || true)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$seq" "$cfg" "$rep" "$jid" "$wall" "$log" "${jrt:-}" >> "$OUT"
}

echo "=== T9 性能对照开始 $(date '+%F %T') ==="
seq=0
for cfg in "${ORDER[@]}"; do
    seq=$((seq + 1))
    run_one "$cfg" "$seq"
done
echo "=== T9 完成 $(date '+%F %T') ==="
echo
echo "=== 结果表 ==="
cat "$OUT"
echo
echo "=== 各配置 pipeline 总时长（跨 rank min/max）==="
for cfg in A B; do
    echo "--- 配置 $cfg ---"
    grep -h "total - complete three-stage pipeline" "$LOGDIR"/${cfg}_${TAG}_PE.run*.log 2>/dev/null \
      | awk '{print $NF, $(NF-1)}' | sort -n | awk 'NR==1{print "  min =",$2,"s"} END{print "  max =",$2,"s"}'
    grep -h "stage 2 - align reads" "$LOGDIR"/${cfg}_${TAG}_PE.run*.log 2>/dev/null \
      | awk '{print $(NF-1)}' | sort -n | awk 'NR==1{print "  stage2 min =",$1,"s"} END{print "  stage2 max =",$1,"s"}'
done
