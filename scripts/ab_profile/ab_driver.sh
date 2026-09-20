#!/usr/bin/env bash
# Drive the A/B + profile runs from the Sunway login node. Each run is its own
# blocking `bsub -I -b` whose direct command is ./SWBWA (see ab_cond.sh for why the
# binary must be the direct child). At most one node is ever in use.
#
#   ./ab_driver.sh smoke    # one cheap condition, validates the whole chain
#   ./ab_driver.sh          # all 12 conditions (safe to re-run: skips complete)
#
# Exported to ab_cond.sh: AB_ROOT. Optional overrides: QUEUE, RETRIES,
# RETRY_SLEEP, DATA_ROOT.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
AB_ROOT=${AB_ROOT:?Set AB_ROOT to the experiment directory}
export AB_ROOT
mkdir -p "$AB_ROOT/logs"
exec 9>"$AB_ROOT/driver.lock"
flock -n 9 || { echo "An A/B driver is already active" >&2; exit 1; }

run_cond()
{
    local cfg=$1 ds=$2 mode=$3
    local out="$AB_ROOT/logs/${cfg}_${ds}_${mode}"
    mkdir -p "$out"

    # Resume support: a condition is complete when all four runs exited 0.
    if [[ -f "$out/run.status" ]]; then
        local ok
        ok=$(awk '$2 == 0' "$out/run.status" | wc -l)
        if [[ "$ok" -eq 4 ]]; then
            echo "[AB] SKIP ${cfg}_${ds}_${mode} (already complete)"
            echo "$cfg $ds $mode 0" >> "$AB_ROOT/logs/conditions.status"
            return 0
        fi
    fi

    echo "[AB] CONDITION START ${cfg}_${ds}_${mode} $(date '+%F %T')"
    local rc=0
    bash "$SCRIPT_DIR/ab_cond.sh" "$cfg" "$ds" "$mode" "$out" || rc=$?
    echo "[AB] CONDITION END   ${cfg}_${ds}_${mode} rc=$rc $(date '+%F %T')"
    echo "$cfg $ds $mode $rc" >> "$AB_ROOT/logs/conditions.status"
    return "$rc"
}

if [[ "${1:-}" == smoke ]]; then
    : > "$AB_ROOT/logs/conditions.status"
    run_cond cross ERR1203383 SE
    echo "[AB] SMOKE DONE $(date '+%F %T')"
    exit 0
fi

: > "$AB_ROOT/logs/conditions.status"

# Cheapest conditions first (cross SE), so usable data lands early.
for cfg in cross single; do
    for mode in SE PE; do
        for ds in ERR1203383 small_SRR7963242 SRR2496709; do
            run_cond "$cfg" "$ds" "$mode"
        done
    done
done

echo "[AB] ALL CONDITIONS DONE $(date '+%F %T')"
