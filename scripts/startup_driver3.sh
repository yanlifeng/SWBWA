#!/usr/bin/env bash
# T8 v3 startup overhead measurement.
#
# Changes from v2 (v2 crashed on $(( )) because polled values went into
# arithmetic):
#   1. the progress poll only writes the raw `bjobs -l` text to a file and
#      does no arithmetic on it;
#   2. the script only records raw epochs (submit / done); all subtraction is
#      done locally afterwards;
#   3. the raw `bjobs -l` text contains both "Submitted from host" and
#      "Started on mn", from which the queueing delay can be computed.
set -euo pipefail

W=${WORK_ROOT:?Set WORK_ROOT to the prepared experiment directory}
DATA_ROOT=${DATA_ROOT:?Set DATA_ROOT to the reference and input directory}
NODES=${NODES:-1}
[[ "$NODES" =~ ^[1-9][0-9]*$ ]] || exit 2
REF=$DATA_ROOT/GRCh38.d1.vd1.fa
R1=$DATA_ROOT/bwa_test_big_data/tiny_1.fastq
R2=$DATA_ROOT/bwa_test_big_data/tiny_2.fastq

LOGDIR=$W/logs/startup3
mkdir -p "$LOGDIR"
exec 9>"$W/driver.lock"
flock -n 9 || { echo "An experiment driver is already active" >&2; exit 1; }
OUT=$LOGDIR/raw_epochs.tsv
[[ ! -e "$OUT" ]] || { echo "Preserving existing results: $OUT" >&2; exit 1; }
printf 'seq\tconfig\trep\tjobid\tsubmit_epoch\tdone_epoch\n' > "$OUT"

run_one() {
    local cfg="$1" rep="$2" seq="$3"
    local rdir args jname log hint bjb out
    local submit_epoch done_epoch jid bpid

    jname="startup3_${cfg}_${rep}"
    if [ "$cfg" = "A" ]; then
        rdir="$W/runs/A_single"
        args=(-N "$NODES" -np 1 -cgsp 64 -mpecg 6 -share_size 2000 -xmalloc -cross_size 42000
              -cache_size 128 -priv_size 16)
    else
        rdir="$W/runs/B_single"
        args=(-N "$NODES" -np 6 -cgsp 64 -share_size 12000 -cache_size 128 -priv_size 16)
    fi

    log="$LOGDIR/$jname.log"
    out="$LOGDIR/$jname.out"
    hint="$W/results/$jname.sam"
    bjb="$LOGDIR/$jname.bjobs_l.txt"
    [[ ! -e "$hint" && ! -e "$out" ]] || {
        echo "Preserving existing files for $jname" >&2; return 1;
    }

    submit_epoch=$(date +%s)
    ( cd "$rdir" && exec bsub -I -b -q q_share "${args[@]}" \
        -J "$jname" -o "$log" \
        ./SWBWA mem -v 4 -t 1 -1 -o "$hint" -I 170,80,500,1 \
        "$REF" "$R1" "$R2" ) < /dev/null > "$out" 2>&1 &
    bpid=$!

    # Bypass polling: persist the raw `bjobs -l` text once it contains
    # "Started on".
    while kill -0 "$bpid" 2>/dev/null; do
        if [ ! -s "$bjb" ]; then
            bjobs -l -J "$jname" > "$bjb.tmp" 2>/dev/null || true
            if grep -q "Started on" "$bjb.tmp" 2>/dev/null; then
                mv -f "$bjb.tmp" "$bjb"
            else
                rm -f "$bjb.tmp"
            fi
        fi
        sleep 60
    done

    wait "$bpid" 2>/dev/null
    grep -q 'has been finished' "$out" || return 1
    grep -q '\[main\] CMD:' "$log" "$out" || return 1
    done_epoch=$(date +%s)

    jid=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' "$out" | head -1)

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$seq" "$cfg" "$rep" "${jid:-}" "$submit_epoch" "$done_epoch" >> "$OUT"

    echo "[T8v3] seq=$seq cfg=$cfg rep=$rep jid=${jid:-?} submit=$submit_epoch done=$done_epoch bjobs_l=$([ -s "$bjb" ] && echo YES || echo NO)"
    rm -f "$hint"
}

echo "=== T8 v3 开始 $(date '+%F %T') ==="
seq=0
for rep in 1 2 3; do
    for cfg in A B; do
        seq=$((seq + 1))
        run_one "$cfg" "$rep" "$seq"
    done
done
echo "=== T8 v3 完成 $(date '+%F %T') ==="
echo
echo "=== 原始 epoch ==="
cat "$OUT"
echo
echo "=== 各次 bjobs -l 中的关键行 ==="
for f in "$LOGDIR"/*.bjobs_l.txt; do
    [ -e "$f" ] || continue
    echo "--- $(basename "$f") ---"
    grep -iE "Submitted from host|Started on" "$f"
done
echo
echo "=== 各次 pipeline 时间 ==="
for f in "$LOGDIR"/startup3_*.log; do
    [ -e "$f" ] || continue
    printf '%s  ' "$(basename "$f")"
    grep -m1 "total - complete three-stage pipeline" "$f" | awk '{print $(NF-1)" s"}'
done
