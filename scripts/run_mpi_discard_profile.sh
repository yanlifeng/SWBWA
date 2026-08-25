#!/usr/bin/env bash
set -euo pipefail

# Environment overrides control the dataset subset, repeat count, and MPI
# layout. MPI_RANKS is the total rank count; bsub -np receives ranks per node.
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}/.."

RESULT_ROOT=${RESULT_ROOT:-correctness_results/mpi_discard_stage2_baseline_20260824}
DATA_ROOT=${DATA_ROOT:-../data/bwa_test_big_data}
REF=${REF:-../data/GRCh38.d1.vd1.fa}
QUEUE=${QUEUE:-q_share}
MPI_RANKS=${MPI_RANKS:-6}
MPI_NODES=${MPI_NODES:-1}
MPI_NODELIST=${MPI_NODELIST:-}
BUILD_JOBS=${BUILD_JOBS:-8}
REPEATS=${REPEATS:-3}

is_positive_integer()
{
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

if ! is_positive_integer "$MPI_NODES" ||
   ! is_positive_integer "$MPI_RANKS" ||
   ! is_positive_integer "$BUILD_JOBS" ||
   ! is_positive_integer "$REPEATS"; then
    echo "MPI_NODES, MPI_RANKS, BUILD_JOBS, and REPEATS must be positive integers" >&2
    exit 1
fi
if ((MPI_RANKS % MPI_NODES != 0)); then
    echo "MPI_RANKS must be divisible by MPI_NODES" >&2
    exit 1
fi
MPI_RANKS_PER_NODE=$((MPI_RANKS / MPI_NODES))
if ((MPI_RANKS_PER_NODE > 6)); then
    echo "MPI ranks per node must not exceed 6" >&2
    exit 1
fi

read -r -a DATASET_LIST <<< \
    "${DATASETS:-ERR1203383 small_SRR7963242 SRR2496709}"
read -r -a READ_MODE_LIST <<< "${READ_MODES:-PE SE}"
for read_mode in "${READ_MODE_LIST[@]}"; do
    case "$read_mode" in
        PE|SE) ;;
        *) echo "READ_MODES accepts only PE and SE" >&2; exit 1 ;;
    esac
done

mkdir -p "$RESULT_ROOT"

bsub_node_args=()
if [[ -n "$MPI_NODELIST" ]]; then
    bsub_node_args=(-node "$MPI_NODELIST")
fi

echo "[BUILD] MPI dynamic + discard, exact read index, single-CG/system"
make clean > "$RESULT_ROOT/build.log" 2>&1
make -j "$BUILD_JOBS" \
    EXEC_MODE=single \
    CPE_ALLOCATOR=system \
    HOST_MALLOC_WRAPPER=1 \
    HOST_MALLOC_STATS=0 \
    CPE_PROFILE=0 \
    USE_MPI=1 \
    MPI_INPUT_MODE=dynamic \
    OUTPUT_MODE=discard \
    MPI_EXACT_READ_INDEX=1 \
    2>&1 | tee -a "$RESULT_ROOT/build.log"

for dataset in "${DATASET_LIST[@]}"; do
    r1="${DATA_ROOT}/${dataset}_1.fastq"
    r2="${DATA_ROOT}/${dataset}_2.fastq"
    [[ -f "$r1" ]] || { echo "missing input: $r1" >&2; exit 1; }
    [[ -f "$r2" ]] || { echo "missing input: $r2" >&2; exit 1; }

    for read_mode in "${READ_MODE_LIST[@]}"; do
        case_logs=()
        for ((repeat = 1; repeat <= REPEATS; ++repeat)); do
            log=$(printf '%s/%s_%s.run%02d.log' \
                         "$RESULT_ROOT" "$dataset" "$read_mode" "$repeat")
            case_logs+=("$log")

            if python3 scripts/check_discard_hash.py "$log" \
                    >/dev/null 2>&1; then
                echo "[SKIP] hash-complete: $log"
                continue
            fi

            output_hint=$(printf '%s/%s_%s.run%02d.discard.sam' \
                                 "$RESULT_ROOT" "$dataset" "$read_mode" \
                                 "$repeat")
            job_name=$(printf 'discard_%s_%s_%02d' \
                              "$dataset" "$read_mode" "$repeat")
            command=(
                env SWBWA_DISCARD_HASH=1 ./SWBWA mem -v 4 -t 1 -1
                -o "$output_hint"
                -I 170,80,500,1
                "$REF" "$r1"
            )
            if [[ "$read_mode" == PE ]]; then
                command+=("$r2")
            fi

            echo "[RUN] $dataset $read_mode repeat=$repeat"
            bsub -I -b -q "$QUEUE" -N "$MPI_NODES" \
                -np "$MPI_RANKS_PER_NODE" \
                "${bsub_node_args[@]}" -cgsp 64 \
                -J "$job_name" -share_size 12000 -cache_size 128 \
                -priv_size 16 -o "$log" "${command[@]}"

            python3 scripts/check_discard_hash.py "$log" >/dev/null
            if [[ -e "$output_hint" ]]; then
                echo "discard mode unexpectedly created $output_hint" >&2
                exit 1
            fi
        done

        python3 scripts/check_discard_hash.py --require-matching \
            "${case_logs[@]}" > \
            "$RESULT_ROOT/${dataset}_${read_mode}.hashes.tsv"
    done
done

python3 scripts/check_discard_hash.py "$RESULT_ROOT"/*.run*.log > \
    "$RESULT_ROOT/output_hashes.tsv"
echo "[DONE] $RESULT_ROOT"
