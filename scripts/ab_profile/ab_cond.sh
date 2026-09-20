#!/usr/bin/env bash
# One A/B condition = 4 blocking bsub jobs, each launching SWBWA as the DIRECT
# command. This is mandatory on Sunway: the CPE/mpe launcher only sets up the CPE
# mapping for the direct child process. Wrapping SWBWA in `bash <script>` makes it
# die with SIGSEGV at startup (the baked-in CPE layout addresses do not match),
# even with a matching data.bin in cwd.
#
# Order per condition: old_np -> new_np -> new_p -> old_p
#   np = no LWPF profiling  -> clean stage2 wall-clock for the A/B pair
#   p  = CPE_PROFILE=1      -> LWPF per-kernel cycle breakdown
#
# The queue is often contended, so a submission can fail immediately with
# "No enough compute nodes"; those are retried with backoff. A job that actually
# ran but returned a non-zero status is NOT retried blindly.
set -u

CFG=${1:?cross|single}; DS=${2:?dataset}; MODE=${3:?PE|SE}; OUT=${4:?out dir}
AB_ROOT=${AB_ROOT:?}
QUEUE=${QUEUE:-q_share}
RETRIES=${RETRIES:-40}
RETRY_SLEEP=${RETRY_SLEEP:-60}
# Root that holds GRCh38.d1.vd1.fa and bwa_test_big_data/.
DATA_ROOT=${DATA_ROOT:?Set DATA_ROOT to the reference and input directory}
[[ "$CFG" == cross || "$CFG" == single ]] || exit 2
[[ "$MODE" == PE || "$MODE" == SE ]] || exit 2

if [[ "$CFG" == cross ]]; then
    ALLOC=(-n 1 -cgsp 64 -share_size 2000 -mpecg 6 -xmalloc -cross_size 42000
           -cache_size 128 -priv_size 16)
else
    ALLOC=(-n 1 -cgsp 64 -share_size 12000 -cache_size 128 -priv_size 16)
fi

REF=$DATA_ROOT/GRCh38.d1.vd1.fa
R1=$DATA_ROOT/bwa_test_big_data/${DS}_1.fastq
R2=$DATA_ROOT/bwa_test_big_data/${DS}_2.fastq
mkdir -p "$OUT" "$AB_ROOT/scratch"
[[ -f "$REF" && -f "$R1" ]] || { echo "[AB] missing ref/input"; exit 2; }
[[ "$MODE" == PE && ! -f "$R2" ]] && { echo "[AB] missing R2"; exit 2; }

node_of()
{
    bjobs -a "$1" 2>/dev/null | awk '
        NR==1 { for (i=1; i<=NF; i++) if ($i == "NODELIST") c = i; next }
        { print $c }' | head -1
}

run_direct()
{
    local label=$1 variant=$2 prof=$3
    local rdir="$AB_ROOT/runs/${variant}_${CFG}_p${prof}"
    local sam="$AB_ROOT/scratch/${CFG}_${DS}_${MODE}.${label}.sam"
    local log="$OUT/${label}.run.log"
    local -a cargs=(mem -v 4 -t 1 -1 -o "$sam" -I 170,80,500,1 "$REF" "$R1")

    [[ "$MODE" == PE ]] && cargs+=("$R2")
    [[ -x "$rdir/SWBWA" ]] || { echo "[AB] missing $rdir/SWBWA"; return 3; }
    [[ ! -e "$sam" ]] || { echo "[AB] preserving existing $sam"; return 3; }

    local attempt rc jid node
    for (( attempt = 1; attempt <= RETRIES; attempt++ )); do
        echo "[AB] submit $CFG/$DS/$MODE $label try=$attempt $(date '+%F %T')"
        ( cd "$rdir" && bsub -I -b -q "$QUEUE" "${ALLOC[@]}" \
            -J "ab_${label}_${CFG}_${DS}_${MODE}" -o "$OUT/${label}.joblog" \
            ./SWBWA "${cargs[@]}" ) > "$log" 2>&1
        rc=$?

        if ! grep -q 'Job <[0-9]' "$log" && grep -q "submit failed" "$log"; then
            echo "[AB] submit rejected ($(tr -d '\n' < "$log")) -> retry in ${RETRY_SLEEP}s"
            sleep "$RETRY_SLEEP"
            continue
        fi

        jid=$(grep -o 'Job <[0-9]*>' "$log" | head -1 | tr -dc '0-9')
        # A detached submission must not cause a second job to be launched.
        if [[ -z "$jid" ]] || ! grep -q 'has been finished' "$log"; then
            echo "[AB] completion unknown; inspect bjobs before restarting" >&2
            return 1
        fi
        if ! grep -q '\[main\] CMD:' "$log" "$OUT/${label}.joblog"; then
            rc=1
        fi
        node="-"
        [[ -n "$jid" ]] && node=$(node_of "$jid")
        [[ -n "$node" ]] || node="-"

        echo "[AB] done   $CFG/$DS/$MODE $label rc=$rc node=$node" \
             "sam_bytes=$(stat -c%s "$sam" 2>/dev/null || echo 0) $(date '+%F %T')"
        echo "$label $rc $node" >> "$OUT/run.status"
        if [[ "$rc" -eq 0 ]]; then rm -f "$sam"; fi
        return $rc
    done

    echo "[AB] GAVE UP $CFG/$DS/$MODE $label after $RETRIES rejected submissions"
    echo "$label 237 -" >> "$OUT/run.status"
    return 237
}

: > "$OUT/run.status"
run_direct old_np old 0 || exit $?
run_direct new_np new 0 || exit $?
run_direct new_p  new 1 || exit $?
run_direct old_p  old 1 || exit $?
echo "[AB] condition complete $CFG/$DS/$MODE $(date '+%F %T')"
