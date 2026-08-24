#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
cd "${SCRIPT_DIR}/.."

EXE=${EXE:-./SWBWA}
REF=${REF:-../data/GRCh38.d1.vd1.fa}
DATA=${DATA:-../data}
RESULT_ROOT=${RESULT_ROOT:-correctness_results}
SORTER=${SORTER:-../fucking_sam_sort}
SORTED_MD5=${SORTED_MD5:-scripts/sorted_sam_md5}
SORTED_MD5_CXX=${SORTED_MD5_CXX:-swg++}
QUEUE=${QUEUE:-q_share}
MPI_NODES=${MPI_NODES:-1}
MPI_RANKS=${MPI_RANKS:-6}
BUILD_JOBS=${BUILD_JOBS:-8}
VERBOSE=${VERBOSE:-4}
SINGLE_SHARE_SIZE=${SINGLE_SHARE_SIZE:-12000}
CGS_SHARE_SIZE=${CGS_SHARE_SIZE:-2000}
MPI_SHARE_SIZE=${MPI_SHARE_SIZE:-12000}
SORT_SHARE_SIZE=${SORT_SHARE_SIZE:-14000}
CACHE_SIZE=${CACHE_SIZE:-128}
PRIVATE_SIZE=${PRIVATE_SIZE:-16}
CROSS_SIZE=${CROSS_SIZE:-42000}
CGS_BSUB_ARGS=(
    -n 1 -cgsp 64
    -share_size "$CGS_SHARE_SIZE" -mpecg 6
    -xmalloc -cross_size "$CROSS_SIZE"
    -cache_size "$CACHE_SIZE" -priv_size "$PRIVATE_SIZE"
)

DEFAULT_DATASETS=(SRR7963242 SRR2496709 ERR1203383)
DEFAULT_BIGDATA_DATASETS=(ERR1203383 small_SRR7963242 SRR2496709)
SELECTED_DATASETS=()
MATRIX_FAILURES=()
MATRIX_FAILURE_COUNT=0
CHUNK_BYTES=""
KEEP_INTERMEDIATE=0
BIGDATA_DATA=${BIGDATA_DATA:-../data/bwa_test_big_data}
BIGDATA_RESULT_ROOT=${BIGDATA_RESULT_ROOT:-bigdata_results}
BIGDATA_MPI_RESULT_ROOT=${BIGDATA_MPI_RESULT_ROOT:-bigdata_mpi_results}
BIGDATA_MD5_FILE=${BIGDATA_MD5_FILE:-scripts/bigdata_expected_md5.tsv}

shopt -s nullglob

usage()
{
    cat <<'EOF'
Usage:
  scripts/correctness.sh run [options]
  scripts/correctness.sh verify [options]
  scripts/correctness.sh matrix [options]
  scripts/correctness.sh bigdata [options]
  scripts/correctness.sh bigdata-mpi [options]

Commands:
  run       Build the selected profile, then run its PE and SE cases.
  verify    Normalize generated SAM files and compare them with known MD5s.
  matrix    Build and test all requested execution and MPI I/O combinations.
  bigdata   Run the final has1/no1 PE/SE regression on three large datasets.
  bigdata-mpi
            Run the four MPI input/output combinations on the large datasets.

Common option:
  --datasets LIST       Comma-separated datasets, or "all" (default: all)

run/matrix options:
  --io-mode MODE        has1, no1, or both (default: has1)
  --chunk-bytes INT     Pass -K INT; omitted by default

run options:
  --profile PROFILE     single, cgs, cgs_cross, or mpi
  --mpi-input MODE      static or dynamic for the MPI profile (default: dynamic)
  --output-mode MODE    ordered, split, or single_unordered
  --output-dir DIR      Directory for generated SAM and log files

verify options:
  --mpi-input MODE      none, static, or dynamic
  --output-mode MODE    ordered, split, or single_unordered
  --output-dir DIR      Directory containing generated SAM files
  --ranks INT           Number of split-output rank files (default: 6)
  --keep-intermediate   Keep concatenated files when verification fails

matrix options:
  --result-root DIR     Root directory for all matrix results
  --keep-intermediate   Keep concatenated files when verification fails

bigdata options:
  --datasets LIST       ERR1203383, small_SRR7963242, SRR2496709, or all
  --result-root DIR     Output root for logs and MD5 summaries
  --chunk-bytes INT     Pass -K INT; omitted by default

bigdata-mpi options:
  --datasets LIST       ERR1203383, small_SRR7963242, SRR2496709, or all
  --result-root DIR     Output root for logs and MD5 summaries
  --chunk-bytes INT     Pass -K INT; omitted by default

The bigdata command always runs has1 and no1 PE/SE tests for:
  single + system
  cgs + system
  cgs + pool
  cgs_cross + system
  cgs_cross + pool
It compares header-free SAM bytes with scripts/bigdata_expected_md5.tsv.
Passing SAM files are removed; failed SAM files are retained for diagnosis.

The bigdata-mpi command runs has1 and no1 for:
  MPI exact index + static  + split
  MPI exact index + static  + single_unordered
  MPI exact index + dynamic + split
  MPI exact index + dynamic + single_unordered
Split files are concatenated in rank order for static input. Dynamic split and
both single_unordered modes are sorted before comparison. Merge and sort work
is submitted to a compute node; at most one job is active at a time.

Passed cases automatically remove all SAM files; logs and MD5 summaries remain.
Matrix runs resume per dataset/read mode from .status checkpoints. Remove the
corresponding result directory when a completely fresh run is required.

Matrix configurations:
  single + system allocator
  cgs_cross + pool allocator
  MPI exact index + static  + split
  MPI exact index + static  + single_unordered
  MPI exact index + dynamic + split
  MPI exact index + dynamic + single_unordered

Environment overrides:
  EXE, REF, DATA, SORTER, QUEUE, MPI_NODES, MPI_RANKS, BUILD_JOBS,
  VERBOSE, SINGLE_SHARE_SIZE, CGS_SHARE_SIZE, MPI_SHARE_SIZE,
  SORT_SHARE_SIZE, CACHE_SIZE, PRIVATE_SIZE, CROSS_SIZE, BIGDATA_DATA,
  BIGDATA_RESULT_ROOT, BIGDATA_MPI_RESULT_ROOT, BIGDATA_MD5_FILE
  SORTED_MD5, SORTED_MD5_CXX

Examples:
  scripts/correctness.sh run --profile single --io-mode has1 \
      --output-dir correctness_results/manual_single
  scripts/correctness.sh verify --mpi-input dynamic --output-mode split \
      --output-dir correctness_results/mpi_dynamic_split/has1
  scripts/correctness.sh matrix --io-mode both
  scripts/correctness.sh bigdata
  scripts/correctness.sh bigdata-mpi
EOF
}

die()
{
    echo "ERROR: $*" >&2
    exit 1
}

validate_positive_integer()
{
    local name=$1
    local value=$2

    [[ "$value" =~ ^[1-9][0-9]*$ ]] ||
        die "$name must be a positive integer: '$value'"
}

normalize_output_mode()
{
    case "$1" in
        ordered|split|single_unordered)
            printf '%s\n' "$1"
            ;;
        single)
            printf '%s\n' single_unordered
            ;;
        *)
            die "invalid output mode '$1'"
            ;;
    esac
}

select_datasets()
{
    local requested=${1:-all}
    local dataset

    SELECTED_DATASETS=()
    if [[ "$requested" == all ]]; then
        SELECTED_DATASETS=("${DEFAULT_DATASETS[@]}")
    else
        IFS=',' read -r -a SELECTED_DATASETS <<< "$requested"
    fi
    ((${#SELECTED_DATASETS[@]} > 0)) || die "no datasets selected"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        case "$dataset" in
            SRR7963242|SRR2496709|ERR1203383)
                ;;
            *)
                die "unknown dataset '$dataset'"
                ;;
        esac
    done
}

select_bigdata_datasets()
{
    local requested=${1:-all}
    local dataset

    SELECTED_DATASETS=()
    if [[ "$requested" == all ]]; then
        SELECTED_DATASETS=("${DEFAULT_BIGDATA_DATASETS[@]}")
    else
        IFS=',' read -r -a SELECTED_DATASETS <<< "$requested"
    fi
    ((${#SELECTED_DATASETS[@]} > 0)) || die "no big-data datasets selected"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        case "$dataset" in
            ERR1203383|small_SRR7963242|SRR2496709)
                ;;
            *)
                die "unknown big-data dataset '$dataset'"
                ;;
        esac
    done
}

select_io_modes()
{
    case "$1" in
        has1)
            IO_MODES=(has1)
            ;;
        no1)
            IO_MODES=(no1)
            ;;
        both)
            IO_MODES=(has1 no1)
            ;;
        *)
            die "invalid I/O mode '$1'; expected has1, no1, or both"
            ;;
    esac
}

expected_md5()
{
    case "$1:$2" in
        SRR7963242:PE) printf '%s\n' a799dc7268f389120ec92820e04b5118 ;;
        SRR7963242:SE) printf '%s\n' 0acfb1f46abd9fed3862c28a33e26da4 ;;
        SRR2496709:PE) printf '%s\n' 20f980ce3955d09e7c133ba26ddf2b77 ;;
        SRR2496709:SE) printf '%s\n' 2ca59c689707cc0aa6e0ce1fe30ae145 ;;
        ERR1203383:PE) printf '%s\n' c2af4bf0b057d5125ce2a9d770f13741 ;;
        ERR1203383:SE) printf '%s\n' 473eec3972fbd651d8b911b9fb5c6e25 ;;
        *) return 1 ;;
    esac
}

bigdata_expected_md5()
{
    local dataset=$1
    local read_mode=$2

    [[ -f "$BIGDATA_MD5_FILE" ]] ||
        die "big-data MD5 reference not found: $BIGDATA_MD5_FILE"
    awk -F '\t' -v dataset="$dataset" -v read_mode="$read_mode" '
        $1 == dataset && $2 == read_mode {
            print $3
            found = 1
            exit
        }
        END {
            if (!found) exit 1
        }
    ' "$BIGDATA_MD5_FILE"
}

file_md5()
{
    local path=$1

    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$path" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then
        md5 -q "$path"
    else
        die "neither md5sum nor md5 is available"
    fi
}

case_paths()
{
    local output_dir=$1
    local dataset=$2
    local read_mode=$3

    CASE_NAME="${dataset}_${read_mode}"
    CASE_BASE="${output_dir}/${CASE_NAME}.sam"
    CASE_LOG="${output_dir}/${CASE_NAME}.run.log"
    CASE_COMBINED="${output_dir}/combined_${CASE_NAME}.sam"
    CASE_NORMALIZED="${output_dir}/normalized_${CASE_NAME}.sam"
    CASE_SORT_LOG="${output_dir}/${CASE_NAME}.sort.log"
    CASE_STATUS_DIR="${output_dir}/.status"
    CASE_RUN_MARKER="${CASE_STATUS_DIR}/${CASE_NAME}.run.ok"
    CASE_MD5_MARKER="${CASE_STATUS_DIR}/${CASE_NAME}.md5.ok"
}

clear_case_outputs()
{
    local output_dir=$1
    local dataset=$2
    local read_mode=$3
    local rank_file

    case_paths "$output_dir" "$dataset" "$read_mode"
    rm -f "$CASE_BASE" "$CASE_LOG" "$CASE_COMBINED" \
          "$CASE_NORMALIZED" "$CASE_SORT_LOG" \
          "$CASE_RUN_MARKER" "$CASE_MD5_MARKER"
    for rank_file in "${CASE_BASE%.sam}".rank*.sam; do
        rm -f "$rank_file"
    done
}

remove_case_sam_outputs()
{
    local output_dir=$1
    local dataset=$2
    local read_mode=$3
    local rank_file

    case_paths "$output_dir" "$dataset" "$read_mode"
    rm -f "$CASE_BASE" "$CASE_COMBINED" "$CASE_NORMALIZED" \
          "$CASE_RUN_MARKER"
    for rank_file in "${CASE_BASE%.sam}".rank*.sam; do
        rm -f "$rank_file"
    done
}

case_outputs_exist()
{
    local output_mode=$1
    local output_dir=$2
    local dataset=$3
    local read_mode=$4
    local ranks=$5
    local rank
    local rank_file

    case_paths "$output_dir" "$dataset" "$read_mode"
    case "$output_mode" in
        ordered|single_unordered)
            [[ -f "$CASE_BASE" ]]
            ;;
        split)
            for ((rank = 0; rank < ranks; ++rank)); do
                rank_file=$(printf '%s.rank%06d.sam' \
                                   "${CASE_BASE%.sam}" "$rank")
                [[ -f "$rank_file" ]] || return 1
            done
            ;;
        *)
            return 1
            ;;
    esac
}

write_run_checkpoint()
{
    local output_mode=$1
    local output_dir=$2
    local dataset=$3
    local read_mode=$4
    local ranks=$5
    local temporary

    case_paths "$output_dir" "$dataset" "$read_mode"
    mkdir -p "$CASE_STATUS_DIR"
    temporary="${CASE_RUN_MARKER}.tmp.$$"
    printf '%s\t%s\n' "$output_mode" "$ranks" > "$temporary"
    mv "$temporary" "$CASE_RUN_MARKER"
}

run_checkpoint_is_valid()
{
    local output_mode=$1
    local output_dir=$2
    local dataset=$3
    local read_mode=$4
    local ranks=$5
    local saved_mode
    local saved_ranks

    case_paths "$output_dir" "$dataset" "$read_mode"
    [[ -f "$CASE_RUN_MARKER" ]] || return 1
    read -r saved_mode saved_ranks < "$CASE_RUN_MARKER" || return 1
    [[ "$saved_mode" == "$output_mode" && "$saved_ranks" == "$ranks" ]] ||
        return 1
    case_outputs_exist "$output_mode" "$output_dir" "$dataset" \
                       "$read_mode" "$ranks"
}

write_md5_checkpoint()
{
    local output_dir=$1
    local dataset=$2
    local read_mode=$3
    local md5=$4
    local temporary

    case_paths "$output_dir" "$dataset" "$read_mode"
    mkdir -p "$CASE_STATUS_DIR"
    temporary="${CASE_MD5_MARKER}.tmp.$$"
    printf '%s\n' "$md5" > "$temporary"
    mv "$temporary" "$CASE_MD5_MARKER"
}

md5_checkpoint_is_valid()
{
    local output_dir=$1
    local dataset=$2
    local read_mode=$3
    local expected
    local actual

    case_paths "$output_dir" "$dataset" "$read_mode"
    [[ -f "$CASE_MD5_MARKER" ]] || return 1
    expected=$(expected_md5 "$dataset" "$read_mode") || return 1
    read -r actual < "$CASE_MD5_MARKER" || return 1
    [[ "$actual" == "$expected" ]]
}

build_mem_command()
{
    local io_mode=$1
    local dataset=$2
    local read_mode=$3
    local output=$4
    local r1="${DATA}/${dataset}_1.fastq_1"
    local r2="${DATA}/${dataset}_2.fastq_1"

    [[ -f "$r1" ]] || die "input file not found: $r1"
    if [[ "$read_mode" == PE ]]; then
        [[ -f "$r2" ]] || die "input file not found: $r2"
    fi
    [[ -f "$REF" ]] || die "reference file not found: $REF"
    [[ -x "$EXE" ]] || die "SWBWA executable not found: $EXE"

    MEM_COMMAND=("$EXE" mem -v "$VERBOSE" -t 1)
    if [[ "$io_mode" == has1 ]]; then
        MEM_COMMAND+=(-1)
    fi
    if [[ -n "$CHUNK_BYTES" ]]; then
        MEM_COMMAND+=(-K "$CHUNK_BYTES")
    fi
    MEM_COMMAND+=(-o "$output" -I 170,80,500,1 "$REF" "$r1")
    if [[ "$read_mode" == PE ]]; then
        MEM_COMMAND+=("$r2")
    fi
}

run_case()
{
    local profile=$1
    local output_mode=$2
    local io_mode=$3
    local output_dir=$4
    local dataset=$5
    local read_mode=$6
    local job_name
    local rank

    mkdir -p "$output_dir"
    clear_case_outputs "$output_dir" "$dataset" "$read_mode"
    case_paths "$output_dir" "$dataset" "$read_mode"
    build_mem_command "$io_mode" "$dataset" "$read_mode" "$CASE_BASE"
    job_name="${profile}_${dataset}_${read_mode}_${io_mode}"

    echo
    echo "============================================================"
    echo "[RUN ] $dataset $read_mode"
    echo "[EXEC] $profile"
    echo "[I/O ] $io_mode"
    echo "[OUT ] $output_mode"
    echo "[SAM ] $CASE_BASE"
    echo "[LOG ] $CASE_LOG"
    echo "============================================================"

    case "$profile" in
        single)
            bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 \
                -J "$job_name" -share_size "$SINGLE_SHARE_SIZE" \
                -o "$CASE_LOG" -cache_size "$CACHE_SIZE" \
                -priv_size "$PRIVATE_SIZE" \
                "${MEM_COMMAND[@]}"
            ;;
        cgs)
            bsub -I -b -q "$QUEUE" "${CGS_BSUB_ARGS[@]}" \
                -J "$job_name" -o "$CASE_LOG" \
                "${MEM_COMMAND[@]}"
            ;;
        cgs_cross)
            bsub -I -b -q "$QUEUE" "${CGS_BSUB_ARGS[@]}" \
                -J "$job_name" -o "$CASE_LOG" \
                "${MEM_COMMAND[@]}"
            ;;
        mpi)
            bsub -I -b -q "$QUEUE" -N "$MPI_NODES" -np "$MPI_RANKS" \
                -cgsp 64 -J "$job_name" -share_size "$MPI_SHARE_SIZE" \
                -o "$CASE_LOG" -cache_size "$CACHE_SIZE" \
                -priv_size "$PRIVATE_SIZE" \
                "${MEM_COMMAND[@]}"
            ;;
        *)
            die "invalid run profile '$profile'"
            ;;
    esac

    if [[ "$output_mode" == split ]]; then
        for ((rank = 0; rank < MPI_RANKS; ++rank)); do
            local rank_file
            rank_file=$(printf '%s.rank%06d.sam' "${CASE_BASE%.sam}" "$rank")
            [[ -f "$rank_file" ]] || die "rank output not generated: $rank_file"
        done
    else
        [[ -f "$CASE_BASE" ]] || die "SAM output not generated: $CASE_BASE"
    fi
    write_run_checkpoint "$output_mode" "$output_dir" "$dataset" \
                         "$read_mode" "$MPI_RANKS"
}

run_selected_cases()
{
    local profile=$1
    local output_mode=$2
    local io_mode=$3
    local output_dir=$4
    local dataset
    local read_mode

    for dataset in "${SELECTED_DATASETS[@]}"; do
        for read_mode in PE SE; do
            run_case "$profile" "$output_mode" "$io_mode" \
                     "$output_dir" "$dataset" "$read_mode"
        done
    done
}

concatenate_rank_outputs()
{
    local base=$1
    local destination=$2
    local ranks=$3
    local temporary="${destination}.tmp.$$"
    local rank
    local rank_file
    local rank_files=()

    for ((rank = 0; rank < ranks; ++rank)); do
        rank_file=$(printf '%s.rank%06d.sam' "${base%.sam}" "$rank")
        [[ -f "$rank_file" ]] || die "missing rank output: $rank_file"
        rank_files+=("$rank_file")
    done

    rm -f "$temporary"
    cat "${rank_files[@]}" > "$temporary"
    mv "$temporary" "$destination"
}

sort_sam()
{
    local input=$1
    local output=$2
    local log=$3
    local job_name=$4

    [[ -x "$SORTER" ]] || die "SAM sorter not found or not executable: $SORTER"
    [[ -f "$input" ]] || die "SAM input not found: $input"
    rm -f "$output" "$log"

    bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 \
        -J "$job_name" -share_size "$SORT_SHARE_SIZE" \
        -o "$log" -cache_size "$CACHE_SIZE" \
        -priv_size "$PRIVATE_SIZE" \
        "$SORTER" "$input" "$output"

    [[ -s "$output" ]] || die "sorted SAM was not generated: $output"
}

verify_case()
{
    local mpi_input=$1
    local output_mode=$2
    local output_dir=$3
    local ranks=$4
    local dataset=$5
    local read_mode=$6
    local source
    local expected
    local actual

    case_paths "$output_dir" "$dataset" "$read_mode"
    expected=$(expected_md5 "$dataset" "$read_mode") ||
        die "no expected MD5 for $dataset $read_mode"

    case "$output_mode" in
        ordered)
            [[ "$mpi_input" == none ]] ||
                die "ordered output is only valid for non-MPI verification"
            [[ -f "$CASE_BASE" ]] || die "SAM output not found: $CASE_BASE"
            source=$CASE_BASE
            ;;
        split)
            [[ "$mpi_input" == static || "$mpi_input" == dynamic ]] ||
                die "split output requires static or dynamic MPI input"
            if [[ "$mpi_input" == static ]]; then
                concatenate_rank_outputs "$CASE_BASE" "$CASE_NORMALIZED" "$ranks"
            else
                concatenate_rank_outputs "$CASE_BASE" "$CASE_COMBINED" "$ranks"
                sort_sam "$CASE_COMBINED" "$CASE_NORMALIZED" \
                         "$CASE_SORT_LOG" "sort_${CASE_NAME}"
                if ((KEEP_INTERMEDIATE == 0)); then
                    rm -f "$CASE_COMBINED"
                fi
            fi
            source=$CASE_NORMALIZED
            ;;
        single_unordered)
            [[ "$mpi_input" == static || "$mpi_input" == dynamic ]] ||
                die "single_unordered output requires MPI input"
            sort_sam "$CASE_BASE" "$CASE_NORMALIZED" \
                     "$CASE_SORT_LOG" "sort_${CASE_NAME}"
            source=$CASE_NORMALIZED
            ;;
        *)
            die "invalid verification output mode '$output_mode'"
            ;;
    esac

    actual=$(file_md5 "$source")
    VERIFY_EXPECTED=$expected
    VERIFY_ACTUAL=$actual
    VERIFY_SOURCE=$source
    if [[ "$actual" == "$expected" ]]; then
        VERIFY_STATUS=PASS
        printf '[PASS] %-12s %-2s  %s\n' "$dataset" "$read_mode" "$actual"
        return 0
    fi

    VERIFY_STATUS=FAIL
    printf '[FAIL] %-12s %-2s  expected=%s actual=%s file=%s\n' \
           "$dataset" "$read_mode" "$expected" "$actual" "$source"
    return 1
}

verify_case_with_checkpoint()
{
    local mpi_input=$1
    local output_mode=$2
    local output_dir=$3
    local ranks=$4
    local dataset=$5
    local read_mode=$6
    local expected
    local actual

    if md5_checkpoint_is_valid "$output_dir" "$dataset" "$read_mode"; then
        case_paths "$output_dir" "$dataset" "$read_mode"
        expected=$(expected_md5 "$dataset" "$read_mode")
        read -r actual < "$CASE_MD5_MARKER"
        VERIFY_STATUS=PASS
        VERIFY_EXPECTED=$expected
        VERIFY_ACTUAL=$actual
        VERIFY_SOURCE="checkpoint:$CASE_MD5_MARKER"
        printf '[SKIP] %-12s %-2s  MD5 already verified (%s)\n' \
               "$dataset" "$read_mode" "$actual"
        return 0
    fi

    if verify_case "$mpi_input" "$output_mode" "$output_dir" "$ranks" \
                   "$dataset" "$read_mode"; then
        remove_case_sam_outputs "$output_dir" "$dataset" "$read_mode"
        write_md5_checkpoint "$output_dir" "$dataset" "$read_mode" \
                              "$VERIFY_ACTUAL"
        printf '[CLEAN] %-12s %-2s  validated SAM files removed\n' \
               "$dataset" "$read_mode"
        return 0
    fi

    case_paths "$output_dir" "$dataset" "$read_mode"
    rm -f "$CASE_RUN_MARKER"
    return 1
}

verify_selected_cases()
{
    local mpi_input=$1
    local output_mode=$2
    local output_dir=$3
    local ranks=$4
    local summary="${output_dir}/md5_results.tsv"
    local failures=0
    local dataset
    local read_mode

    mkdir -p "$output_dir"
    printf 'dataset\tread_mode\tstatus\texpected\tactual\tchecked_file\n' > "$summary"

    echo
    echo "============================================================"
    echo "[VERIFY] mpi_input=$mpi_input output_mode=$output_mode"
    echo "[DIR   ] $output_dir"
    echo "============================================================"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        for read_mode in PE SE; do
            if verify_case_with_checkpoint "$mpi_input" "$output_mode" \
                                           "$output_dir" "$ranks" \
                                           "$dataset" "$read_mode"; then
                :
            else
                failures=$((failures + 1))
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                   "$dataset" "$read_mode" "$VERIFY_STATUS" \
                   "$VERIFY_EXPECTED" "$VERIFY_ACTUAL" "$VERIFY_SOURCE" \
                   >> "$summary"
        done
    done

    echo "[SUMMARY] $summary"
    ((failures == 0))
}

run_and_verify_selected_cases()
{
    local profile=$1
    local mpi_input=$2
    local output_mode=$3
    local io_mode=$4
    local output_dir=$5
    local ranks=$6
    local summary="${output_dir}/md5_results.tsv"
    local failures=0
    local dataset
    local read_mode

    mkdir -p "$output_dir"
    printf 'dataset\tread_mode\tstatus\texpected\tactual\tchecked_file\n' \
        > "$summary"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        for read_mode in PE SE; do
            if ! md5_checkpoint_is_valid "$output_dir" "$dataset" \
                                           "$read_mode"; then
                if run_checkpoint_is_valid "$output_mode" "$output_dir" \
                                           "$dataset" "$read_mode" "$ranks"; then
                    printf '[RESUME] %-12s %-2s  alignment complete; continuing verification\n' \
                           "$dataset" "$read_mode"
                else
                    run_case "$profile" "$output_mode" "$io_mode" \
                             "$output_dir" "$dataset" "$read_mode"
                fi
            fi

            if verify_case_with_checkpoint "$mpi_input" "$output_mode" \
                                           "$output_dir" "$ranks" \
                                           "$dataset" "$read_mode"; then
                :
            else
                failures=$((failures + 1))
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                   "$dataset" "$read_mode" "$VERIFY_STATUS" \
                   "$VERIFY_EXPECTED" "$VERIFY_ACTUAL" "$VERIFY_SOURCE" \
                   >> "$summary"
        done
    done

    echo "[SUMMARY] $summary"
    ((failures == 0))
}

configuration_requires_alignment()
{
    local config_dir=$1
    local output_mode=$2
    local io_mode
    local output_dir
    local dataset
    local read_mode

    for io_mode in "${IO_MODES[@]}"; do
        output_dir="${config_dir}/${io_mode}"
        for dataset in "${SELECTED_DATASETS[@]}"; do
            for read_mode in PE SE; do
                if md5_checkpoint_is_valid "$output_dir" "$dataset" \
                                            "$read_mode"; then
                    continue
                fi
                if ! run_checkpoint_is_valid "$output_mode" "$output_dir" \
                                             "$dataset" "$read_mode" \
                                             "$MPI_RANKS"; then
                    return 0
                fi
            done
        done
    done
    return 1
}

build_configuration()
{
    local output_dir=$1
    local exec_mode=$2
    local allocator=$3
    local use_mpi=$4
    local mpi_input=$5
    local output_mode=$6
    local build_log="${output_dir}/build.log"
    local make_args=(
        EXEC_MODE="$exec_mode"
        CPE_ALLOCATOR="$allocator"
        HOST_MALLOC_WRAPPER=1
        HOST_MALLOC_STATS=0
        CPE_PROFILE=0
        USE_MPI="$use_mpi"
    )

    if [[ "$use_mpi" == 1 ]]; then
        make_args+=(
            MPI_INPUT_MODE="$mpi_input"
            OUTPUT_MODE="$output_mode"
            MPI_EXACT_READ_INDEX=1
        )
    fi

    mkdir -p "$output_dir"
    echo "[BUILD] EXEC_MODE=$exec_mode CPE_ALLOCATOR=$allocator USE_MPI=$use_mpi"
    if [[ "$use_mpi" == 1 ]]; then
        echo "[BUILD] MPI_INPUT_MODE=$mpi_input OUTPUT_MODE=$output_mode EXACT=1"
    fi
    if [[ "$exec_mode" == cgs_cross ]]; then
        [[ -x ./build_cross.sh ]] || die "build_cross.sh is not executable"
        echo "[BUILD] DRIVER=build_cross.sh (two-pass CPE layout build)"
        ./build_cross.sh "$BUILD_JOBS" "${make_args[@]}" \
            2>&1 | tee "$build_log"
    else
        echo "[BUILD] DRIVER=make"
        make clean > "$build_log" 2>&1
        make -j "$BUILD_JOBS" "${make_args[@]}" \
            2>&1 | tee -a "$build_log"
    fi
}

run_matrix_configuration()
{
    local name=$1
    local profile=$2
    local exec_mode=$3
    local allocator=$4
    local use_mpi=$5
    local mpi_input=$6
    local output_mode=$7
    local config_dir="${RESULT_ROOT}/${name}"
    local io_mode
    local run_dir

    echo
    echo "################################################################"
    echo "# MATRIX CONFIGURATION: $name"
    echo "################################################################"
    if configuration_requires_alignment "$config_dir" "$output_mode"; then
        build_configuration "$config_dir" "$exec_mode" "$allocator" \
                            "$use_mpi" "$mpi_input" "$output_mode"
    else
        echo "[BUILD SKIP] all remaining cases only need verification or are complete"
    fi

    for io_mode in "${IO_MODES[@]}"; do
        run_dir="${config_dir}/${io_mode}"
        if run_and_verify_selected_cases "$profile" "$mpi_input" \
                                         "$output_mode" "$io_mode" \
                                         "$run_dir" "$MPI_RANKS"; then
            echo "[CONFIG PASS] $name $io_mode"
        else
            echo "[CONFIG FAIL] $name $io_mode"
            MATRIX_FAILURES+=("${name}:${io_mode}")
            MATRIX_FAILURE_COUNT=$((MATRIX_FAILURE_COUNT + 1))
        fi
    done
}

run_command()
{
    local profile=""
    local io_mode=has1
    local mpi_input=""
    local output_mode=""
    local output_dir=""
    local datasets=all
    local exec_mode
    local allocator
    local use_mpi
    local io
    local io_output_dir

    while (($# > 0)); do
        case "$1" in
            --profile) profile=${2:?}; shift 2 ;;
            --io-mode) io_mode=${2:?}; shift 2 ;;
            --mpi-input) mpi_input=${2:?}; shift 2 ;;
            --output-mode) output_mode=${2:?}; shift 2 ;;
            --output-dir) output_dir=${2:?}; shift 2 ;;
            --datasets) datasets=${2:?}; shift 2 ;;
            --chunk-bytes) CHUNK_BYTES=${2:?}; shift 2 ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown run option '$1'" ;;
        esac
    done

    [[ -n "$profile" ]] || die "run requires --profile"
    case "$profile" in
        single)
            exec_mode=single
            allocator=system
            use_mpi=0
            output_mode=${output_mode:-ordered}
            ;;
        cgs)
            exec_mode=cgs
            allocator=system
            use_mpi=0
            output_mode=${output_mode:-ordered}
            ;;
        cgs_cross)
            exec_mode=cgs_cross
            allocator=pool
            use_mpi=0
            output_mode=${output_mode:-ordered}
            ;;
        mpi)
            exec_mode=single
            allocator=system
            use_mpi=1
            mpi_input=${mpi_input:-dynamic}
            output_mode=${output_mode:-single_unordered}
            ;;
        *)
            die "invalid profile '$profile'"
            ;;
    esac
    output_mode=$(normalize_output_mode "$output_mode")
    if [[ "$profile" != mpi && "$output_mode" != ordered ]]; then
        die "non-MPI profiles require --output-mode ordered"
    fi
    if [[ "$profile" != mpi && -n "$mpi_input" ]]; then
        die "--mpi-input is only valid with --profile mpi"
    fi
    if [[ "$profile" == mpi && "$output_mode" == ordered ]]; then
        die "MPI profile requires split or single_unordered output"
    fi
    if [[ "$profile" == mpi ]]; then
        case "$mpi_input" in
            static|dynamic) ;;
            *) die "invalid MPI input '$mpi_input'" ;;
        esac
    fi
    if [[ -n "$CHUNK_BYTES" ]]; then
        validate_positive_integer --chunk-bytes "$CHUNK_BYTES"
    fi
    validate_positive_integer MPI_RANKS "$MPI_RANKS"
    select_datasets "$datasets"
    select_io_modes "$io_mode"
    output_dir=${output_dir:-${RESULT_ROOT}/manual_${profile}}

    build_configuration "$output_dir" "$exec_mode" "$allocator" \
                        "$use_mpi" "$mpi_input" "$output_mode"

    for io in "${IO_MODES[@]}"; do
        io_output_dir=$output_dir
        if ((${#IO_MODES[@]} > 1)); then
            io_output_dir="${output_dir}/${io}"
        fi
        run_selected_cases "$profile" "$output_mode" "$io" "$io_output_dir"
    done
}

verify_command()
{
    local mpi_input=none
    local output_mode=ordered
    local output_dir=""
    local datasets=all
    local ranks=$MPI_RANKS

    while (($# > 0)); do
        case "$1" in
            --mpi-input) mpi_input=${2:?}; shift 2 ;;
            --output-mode) output_mode=${2:?}; shift 2 ;;
            --output-dir) output_dir=${2:?}; shift 2 ;;
            --datasets) datasets=${2:?}; shift 2 ;;
            --ranks) ranks=${2:?}; shift 2 ;;
            --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown verify option '$1'" ;;
        esac
    done

    [[ -n "$output_dir" ]] || die "verify requires --output-dir"
    case "$mpi_input" in none|static|dynamic) ;; *) die "invalid MPI input '$mpi_input'" ;; esac
    output_mode=$(normalize_output_mode "$output_mode")
    validate_positive_integer --ranks "$ranks"
    select_datasets "$datasets"
    verify_selected_cases "$mpi_input" "$output_mode" "$output_dir" "$ranks"
}

matrix_command()
{
    local io_mode=has1
    local datasets=all

    while (($# > 0)); do
        case "$1" in
            --io-mode) io_mode=${2:?}; shift 2 ;;
            --datasets) datasets=${2:?}; shift 2 ;;
            --result-root) RESULT_ROOT=${2:?}; shift 2 ;;
            --chunk-bytes) CHUNK_BYTES=${2:?}; shift 2 ;;
            --keep-intermediate) KEEP_INTERMEDIATE=1; shift ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown matrix option '$1'" ;;
        esac
    done

    if [[ -n "$CHUNK_BYTES" ]]; then
        validate_positive_integer --chunk-bytes "$CHUNK_BYTES"
    fi
    validate_positive_integer MPI_RANKS "$MPI_RANKS"
    validate_positive_integer BUILD_JOBS "$BUILD_JOBS"
    select_datasets "$datasets"
    select_io_modes "$io_mode"
    EXE=./SWBWA
    mkdir -p "$RESULT_ROOT"

    run_matrix_configuration single_system single single system 0 none ordered
    run_matrix_configuration cgs_cross_pool cgs_cross cgs_cross pool 0 none ordered
    run_matrix_configuration mpi_static_split mpi single system 1 static split
    run_matrix_configuration mpi_static_single mpi single system 1 static single_unordered
    run_matrix_configuration mpi_dynamic_split mpi single system 1 dynamic split
    run_matrix_configuration mpi_dynamic_single mpi single system 1 dynamic single_unordered

    echo
    echo "################################################################"
    echo "# MATRIX SUMMARY"
    echo "################################################################"
    if ((MATRIX_FAILURE_COUNT == 0)); then
        echo "All correctness configurations passed."
        return 0
    fi

    echo "Failed configurations:"
    printf '  %s\n' "${MATRIX_FAILURES[@]}"
    return 1
}

bigdata_case_paths()
{
    local config_dir=$1
    local dataset=$2
    local read_mode=$3

    BIGDATA_CASE_NAME="${dataset}_${read_mode}"
    BIGDATA_CASE_SAM="${config_dir}/${BIGDATA_CASE_NAME}.sam"
    BIGDATA_CASE_LOG="${config_dir}/${BIGDATA_CASE_NAME}.run.log"
    BIGDATA_CASE_COMBINED="${config_dir}/combined_${BIGDATA_CASE_NAME}.sam"
    BIGDATA_CASE_NORMALIZED="${config_dir}/normalized_${BIGDATA_CASE_NAME}.sam"
    BIGDATA_CASE_MERGE_LOG="${config_dir}/${BIGDATA_CASE_NAME}.merge.log"
    BIGDATA_CASE_SORT_LOG="${config_dir}/${BIGDATA_CASE_NAME}.sort.log"
    BIGDATA_CASE_CONCAT_MD5="${config_dir}/${BIGDATA_CASE_NAME}.concat.md5"
    BIGDATA_CASE_SORTED_MD5="${config_dir}/${BIGDATA_CASE_NAME}.sorted.md5"
    BIGDATA_STATUS_DIR="${config_dir}/.status"
    BIGDATA_MD5_MARKER="${BIGDATA_STATUS_DIR}/${BIGDATA_CASE_NAME}.md5.ok"
}

bigdata_checkpoint_is_valid()
{
    local config_dir=$1
    local dataset=$2
    local read_mode=$3
    local expected
    local actual

    bigdata_case_paths "$config_dir" "$dataset" "$read_mode"
    [[ -f "$BIGDATA_MD5_MARKER" ]] || return 1
    expected=$(bigdata_expected_md5 "$dataset" "$read_mode") || return 1
    read -r actual < "$BIGDATA_MD5_MARKER" || return 1
    [[ "$actual" == "$expected" ]]
}

bigdata_configuration_requires_run()
{
    local config_dir=$1
    local io_mode
    local run_dir
    local dataset
    local read_mode

    for io_mode in "${IO_MODES[@]}"; do
        run_dir="${config_dir}/${io_mode}"
        for dataset in "${SELECTED_DATASETS[@]}"; do
            for read_mode in PE SE; do
                if ! bigdata_checkpoint_is_valid "$run_dir" "$dataset" \
                                                    "$read_mode"; then
                    return 0
                fi
            done
        done
    done
    return 1
}

build_bigdata_mem_command()
{
    local io_mode=$1
    local dataset=$2
    local read_mode=$3
    local output=$4
    local r1="${BIGDATA_DATA}/${dataset}_1.fastq"
    local r2="${BIGDATA_DATA}/${dataset}_2.fastq"

    [[ -f "$r1" ]] || die "big-data input file not found: $r1"
    if [[ "$read_mode" == PE ]]; then
        [[ -f "$r2" ]] || die "big-data input file not found: $r2"
    fi
    [[ -f "$REF" ]] || die "reference file not found: $REF"
    [[ -x "$EXE" ]] || die "SWBWA executable not found: $EXE"

    BIGDATA_MEM_COMMAND=("$EXE" mem -v "$VERBOSE" -t 1)
    if [[ "$io_mode" == has1 ]]; then
        BIGDATA_MEM_COMMAND+=(-1)
    fi
    if [[ -n "$CHUNK_BYTES" ]]; then
        BIGDATA_MEM_COMMAND+=(-K "$CHUNK_BYTES")
    fi
    BIGDATA_MEM_COMMAND+=(-o "$output" -I 170,80,500,1 "$REF" "$r1")
    if [[ "$read_mode" == PE ]]; then
        BIGDATA_MEM_COMMAND+=("$r2")
    fi
}

submit_bigdata_job()
{
    local profile=$1
    local job_name=$2
    local log=$3

    case "$profile" in
        single)
            bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 \
                -J "$job_name" -share_size "$SINGLE_SHARE_SIZE" \
                -o "$log" -cache_size "$CACHE_SIZE" \
                -priv_size "$PRIVATE_SIZE" \
                "${BIGDATA_MEM_COMMAND[@]}"
            ;;
        cgs)
            bsub -I -b -q "$QUEUE" "${CGS_BSUB_ARGS[@]}" \
                -J "$job_name" -o "$log" \
                "${BIGDATA_MEM_COMMAND[@]}"
            ;;
        cgs_cross)
            bsub -I -b -q "$QUEUE" "${CGS_BSUB_ARGS[@]}" \
                -J "$job_name" -o "$log" \
                "${BIGDATA_MEM_COMMAND[@]}"
            ;;
        mpi)
            bsub -I -b -q "$QUEUE" -N 1 -np "$MPI_RANKS" \
                -cgsp 64 -J "$job_name" -share_size "$MPI_SHARE_SIZE" \
                -o "$log" -cache_size "$CACHE_SIZE" \
                -priv_size "$PRIVATE_SIZE" \
                "${BIGDATA_MEM_COMMAND[@]}"
            ;;
        *)
            die "invalid big-data run profile '$profile'"
            ;;
    esac
}

remove_bigdata_case_sams()
{
    local config_dir=$1
    local dataset=$2
    local read_mode=$3
    local rank_file

    bigdata_case_paths "$config_dir" "$dataset" "$read_mode"
    rm -f "$BIGDATA_CASE_SAM" "$BIGDATA_CASE_COMBINED" \
          "$BIGDATA_CASE_NORMALIZED"
    for rank_file in "${BIGDATA_CASE_SAM%.sam}".rank*.sam; do
        rm -f "$rank_file"
    done
}

md5_bigdata_rank_outputs()
{
    local base=$1
    local result=$2
    local log=$3
    local job_name=$4
    local rank
    local rank_file
    local rank_files=()

    for ((rank = 0; rank < MPI_RANKS; ++rank)); do
        rank_file=$(printf '%s.rank%06d.sam' "${base%.sam}" "$rank")
        [[ -f "$rank_file" ]] || die "missing rank output: $rank_file"
        rank_files+=("$rank_file")
    done

    rm -f "$result" "$log"
    if ! bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 \
        -J "$job_name" -share_size "$SORT_SHARE_SIZE" \
        -o "$log" -cache_size "$CACHE_SIZE" \
        -priv_size "$PRIVATE_SIZE" \
        scripts/md5_concat_sam.sh "$result" "${rank_files[@]}"; then
        [[ -s "$result" ]] || return 1
    fi
    [[ -s "$result" ]] || die "concatenated SAM MD5 was not generated: $result"
}

collect_bigdata_rank_outputs()
{
    local base=$1
    local rank
    local rank_file

    BIGDATA_RANK_FILES=()
    for ((rank = 0; rank < MPI_RANKS; ++rank)); do
        rank_file=$(printf '%s.rank%06d.sam' "${base%.sam}" "$rank")
        [[ -f "$rank_file" ]] || die "missing rank output: $rank_file"
        BIGDATA_RANK_FILES+=("$rank_file")
    done
}

sorted_bigdata_md5()
{
    local result=$1
    local log=$2
    local job_name=$3
    shift 3

    [[ -x "$SORTED_MD5" ]] ||
        die "sorted-SAM MD5 tool not found or not executable: $SORTED_MD5"
    (($# > 0)) || die "sorted-SAM MD5 requires at least one SAM input"
    rm -f "$result" "$log"
    bsub -I -b -q "$QUEUE" -n 1 -cgsp 64 \
        -J "$job_name" -share_size "$SORT_SHARE_SIZE" \
        -o "$log" -cache_size "$CACHE_SIZE" \
        -priv_size "$PRIVATE_SIZE" \
        "$SORTED_MD5" "$result" "$@"
    [[ -s "$result" ]] || die "sorted SAM MD5 was not generated: $result"
}

build_sorted_md5_tool()
{
    local source=scripts/sorted_sam_md5.cpp

    [[ -f "$source" ]] || die "sorted-SAM MD5 source not found: $source"
    if [[ ! -x "$SORTED_MD5" || "$source" -nt "$SORTED_MD5" ]]; then
        echo "[BUILD] sorted-SAM MD5 helper: $SORTED_MD5"
        "$SORTED_MD5_CXX" -mhybrid -static -O2 -std=c++11 \
            -D_FILE_OFFSET_BITS=64 \
            "$source" -o "$SORTED_MD5"
    fi
}

read_md5_result()
{
    local path=$1
    local attempt
    local value

    for ((attempt = 0; attempt < 30; ++attempt)); do
        value=$(awk 'NR == 1 { print $1; exit }' "$path" 2>/dev/null || true)
        if [[ "$value" =~ ^[0-9a-fA-F]{32}$ ]]; then
            printf '%s\n' "$value"
            return 0
        fi
        sleep 1
    done
    return 1
}

verify_bigdata_mpi_case()
{
    local config_dir=$1
    local dataset=$2
    local read_mode=$3
    local mpi_input=$4
    local output_mode=$5
    local expected
    local actual
    local source

    bigdata_case_paths "$config_dir" "$dataset" "$read_mode"
    expected=$(bigdata_expected_md5 "$dataset" "$read_mode") ||
        die "no big-data expected MD5 for $dataset $read_mode"

    case "$output_mode" in
        split)
            if [[ "$mpi_input" == static ]]; then
                md5_bigdata_rank_outputs \
                    "$BIGDATA_CASE_SAM" "$BIGDATA_CASE_CONCAT_MD5" \
                    "$BIGDATA_CASE_MERGE_LOG" "md5_${BIGDATA_CASE_NAME}"
                actual=$(read_md5_result "$BIGDATA_CASE_CONCAT_MD5") ||
                    die "invalid concatenated SAM MD5: $BIGDATA_CASE_CONCAT_MD5"
                source="stream:${BIGDATA_CASE_CONCAT_MD5}"
            else
                collect_bigdata_rank_outputs "$BIGDATA_CASE_SAM"
                sorted_bigdata_md5 "$BIGDATA_CASE_SORTED_MD5" \
                    "$BIGDATA_CASE_SORT_LOG" "sort_${BIGDATA_CASE_NAME}" \
                    "${BIGDATA_RANK_FILES[@]}"
                actual=$(read_md5_result "$BIGDATA_CASE_SORTED_MD5") ||
                    die "invalid sorted SAM MD5: $BIGDATA_CASE_SORTED_MD5"
                source="sorted-stream:${BIGDATA_CASE_SORTED_MD5}"
            fi
            ;;
        single_unordered)
            [[ -f "$BIGDATA_CASE_SAM" ]] ||
                die "SAM output not found: $BIGDATA_CASE_SAM"
            sorted_bigdata_md5 "$BIGDATA_CASE_SORTED_MD5" \
                "$BIGDATA_CASE_SORT_LOG" "sort_${BIGDATA_CASE_NAME}" \
                "$BIGDATA_CASE_SAM"
            actual=$(read_md5_result "$BIGDATA_CASE_SORTED_MD5") ||
                die "invalid sorted SAM MD5: $BIGDATA_CASE_SORTED_MD5"
            source="sorted-stream:${BIGDATA_CASE_SORTED_MD5}"
            ;;
        *)
            die "invalid MPI big-data output mode '$output_mode'"
            ;;
    esac

    if [[ -z "${actual:-}" ]]; then
        actual=$(file_md5 "$source")
    fi
    BIGDATA_VERIFY_EXPECTED=$expected
    BIGDATA_VERIFY_ACTUAL=$actual
    if [[ "$actual" == "$expected" ]]; then
        remove_bigdata_case_sams "$config_dir" "$dataset" "$read_mode"
        mkdir -p "$BIGDATA_STATUS_DIR"
        printf '%s\n' "$actual" > "$BIGDATA_MD5_MARKER"
        BIGDATA_VERIFY_STATUS=PASS
        printf '[PASS] %-18s %-2s  %s  (SAM files removed)\n' \
               "$dataset" "$read_mode" "$actual"
        return 0
    fi

    BIGDATA_VERIFY_STATUS=FAIL
    printf '[FAIL] %-18s %-2s  expected=%s actual=%s file=%s\n' \
           "$dataset" "$read_mode" "$expected" "$actual" "$source"
    return 1
}

run_bigdata_mpi_case()
{
    local configuration=$1
    local io_mode=$2
    local run_dir=$3
    local dataset=$4
    local read_mode=$5
    local mpi_input=$6
    local output_mode=$7
    local expected
    local actual="-"
    local rank
    local rank_file
    local job_name="big_${configuration}_${dataset}_${read_mode}_${io_mode}"

    bigdata_case_paths "$run_dir" "$dataset" "$read_mode"
    expected=$(bigdata_expected_md5 "$dataset" "$read_mode") ||
        die "no big-data expected MD5 for $dataset $read_mode"
    remove_bigdata_case_sams "$run_dir" "$dataset" "$read_mode"
    rm -f "$BIGDATA_CASE_LOG" "$BIGDATA_CASE_MERGE_LOG" \
          "$BIGDATA_CASE_SORT_LOG" "$BIGDATA_CASE_CONCAT_MD5" \
          "$BIGDATA_CASE_SORTED_MD5" "$BIGDATA_MD5_MARKER"
    build_bigdata_mem_command "$io_mode" "$dataset" "$read_mode" \
                              "$BIGDATA_CASE_SAM"

    echo
    echo "============================================================"
    echo "[BIG ] $dataset $read_mode"
    echo "[CONF] $configuration"
    echo "[I/O ] $io_mode"
    echo "[MPI ] input=$mpi_input output=$output_mode exact_index=1"
    echo "[SAM ] $BIGDATA_CASE_SAM"
    echo "[LOG ] $BIGDATA_CASE_LOG"
    echo "============================================================"

    if ! submit_bigdata_job mpi "$job_name" "$BIGDATA_CASE_LOG"; then
        BIGDATA_VERIFY_STATUS=ERROR
        BIGDATA_VERIFY_EXPECTED=$expected
        BIGDATA_VERIFY_ACTUAL=$actual
        printf '[ERROR] %-18s %-2s  SWBWA job failed; partial SAM retained\n' \
               "$dataset" "$read_mode"
        return 1
    fi

    if [[ "$output_mode" == split ]]; then
        for ((rank = 0; rank < MPI_RANKS; ++rank)); do
            rank_file=$(printf '%s.rank%06d.sam' \
                               "${BIGDATA_CASE_SAM%.sam}" "$rank")
            if [[ ! -f "$rank_file" ]]; then
                BIGDATA_VERIFY_STATUS=ERROR
                BIGDATA_VERIFY_EXPECTED=$expected
                BIGDATA_VERIFY_ACTUAL=$actual
                printf '[ERROR] missing rank output: %s\n' "$rank_file"
                return 1
            fi
        done
    elif [[ ! -f "$BIGDATA_CASE_SAM" ]]; then
        BIGDATA_VERIFY_STATUS=ERROR
        BIGDATA_VERIFY_EXPECTED=$expected
        BIGDATA_VERIFY_ACTUAL=$actual
        printf '[ERROR] SAM output was not generated: %s\n' "$BIGDATA_CASE_SAM"
        return 1
    fi

    verify_bigdata_mpi_case "$run_dir" "$dataset" "$read_mode" \
                            "$mpi_input" "$output_mode"
}

run_bigdata_mpi_configuration()
{
    local configuration=$1
    local mpi_input=$2
    local output_mode=$3
    local config_dir="${BIGDATA_MPI_RESULT_ROOT}/${configuration}"
    local summary="${config_dir}/md5_results.tsv"
    local failures=0
    local io_mode
    local run_dir
    local dataset
    local read_mode

    echo
    echo "################################################################"
    echo "# BIG-DATA MPI CONFIGURATION: $configuration"
    echo "################################################################"
    mkdir -p "$config_dir"
    if bigdata_configuration_requires_run "$config_dir"; then
        if ! build_configuration "$config_dir" single system 1 \
                                 "$mpi_input" "$output_mode"; then
            echo "[BUILD FAIL] $configuration" >&2
            return 1
        fi
    else
        echo "[BUILD SKIP] all selected big-data MPI cases already passed"
    fi

    printf 'configuration\tio_mode\tmpi_input\toutput_mode\tdataset\tread_mode\tstatus\texpected\tactual\tlog\n' \
        > "$summary"
    for io_mode in "${IO_MODES[@]}"; do
        run_dir="${config_dir}/${io_mode}"
        mkdir -p "$run_dir"
        for dataset in "${SELECTED_DATASETS[@]}"; do
            for read_mode in PE SE; do
                bigdata_case_paths "$run_dir" "$dataset" "$read_mode"
                if bigdata_checkpoint_is_valid "$run_dir" "$dataset" \
                                                   "$read_mode"; then
                    BIGDATA_VERIFY_STATUS=PASS
                    BIGDATA_VERIFY_EXPECTED=$(bigdata_expected_md5 \
                        "$dataset" "$read_mode")
                    read -r BIGDATA_VERIFY_ACTUAL < "$BIGDATA_MD5_MARKER"
                    printf '[SKIP] %-5s %-18s %-2s  MD5 already verified (%s)\n' \
                           "$io_mode" "$dataset" "$read_mode" \
                           "$BIGDATA_VERIFY_ACTUAL"
                elif ! run_bigdata_mpi_case "$configuration" "$io_mode" \
                                              "$run_dir" "$dataset" \
                                              "$read_mode" "$mpi_input" \
                                              "$output_mode"; then
                    failures=$((failures + 1))
                fi
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                       "$configuration" "$io_mode" "$mpi_input" \
                       "$output_mode" "$dataset" "$read_mode" \
                       "$BIGDATA_VERIFY_STATUS" \
                       "$BIGDATA_VERIFY_EXPECTED" \
                       "$BIGDATA_VERIFY_ACTUAL" "$BIGDATA_CASE_LOG" \
                       >> "$summary"
            done
        done
    done
    echo "[SUMMARY] $summary"
    ((failures == 0))
}

write_bigdata_mpi_global_summary()
{
    local output="${BIGDATA_MPI_RESULT_ROOT}/md5_results.tsv"
    local configuration
    local summary

    printf 'configuration\tio_mode\tmpi_input\toutput_mode\tdataset\tread_mode\tstatus\texpected\tactual\tlog\n' \
        > "$output"
    for configuration in mpi_static_split mpi_static_single \
                         mpi_dynamic_split mpi_dynamic_single; do
        summary="${BIGDATA_MPI_RESULT_ROOT}/${configuration}/md5_results.tsv"
        [[ -f "$summary" ]] || continue
        tail -n +2 "$summary" >> "$output"
    done
    echo "[GLOBAL SUMMARY] $output"
}

bigdata_mpi_command()
{
    local datasets=all
    local failures=0

    while (($# > 0)); do
        case "$1" in
            --datasets) datasets=${2:?}; shift 2 ;;
            --result-root) BIGDATA_MPI_RESULT_ROOT=${2:?}; shift 2 ;;
            --chunk-bytes) CHUNK_BYTES=${2:?}; shift 2 ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown bigdata-mpi option '$1'" ;;
        esac
    done

    if [[ -n "$CHUNK_BYTES" ]]; then
        validate_positive_integer --chunk-bytes "$CHUNK_BYTES"
    fi
    validate_positive_integer BUILD_JOBS "$BUILD_JOBS"
    validate_positive_integer MPI_RANKS "$MPI_RANKS"
    [[ "$MPI_NODES" == 1 ]] ||
        die "bigdata-mpi is restricted to MPI_NODES=1"
    select_bigdata_datasets "$datasets"
    select_io_modes both
    [[ -d "$BIGDATA_DATA" ]] ||
        die "big-data input directory not found: $BIGDATA_DATA"
    [[ -f "$BIGDATA_MD5_FILE" ]] ||
        die "big-data MD5 reference not found: $BIGDATA_MD5_FILE"
    EXE=./SWBWA
    mkdir -p "$BIGDATA_MPI_RESULT_ROOT"
    build_sorted_md5_tool

    if ! run_bigdata_mpi_configuration mpi_static_split static split; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_mpi_configuration mpi_static_single static \
                                        single_unordered; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_mpi_configuration mpi_dynamic_split dynamic split; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_mpi_configuration mpi_dynamic_single dynamic \
                                        single_unordered; then
        failures=$((failures + 1))
    fi
    write_bigdata_mpi_global_summary

    echo
    if ((failures == 0)); then
        echo "All big-data MPI correctness configurations passed."
        return 0
    fi
    echo "$failures big-data MPI configuration(s) had failures."
    return 1
}

run_bigdata_case()
{
    local configuration=$1
    local profile=$2
    local io_mode=$3
    local run_dir=$4
    local dataset=$5
    local read_mode=$6
    local expected
    local actual="-"
    local job_name="big_${configuration}_${dataset}_${read_mode}_${io_mode}"

    bigdata_case_paths "$run_dir" "$dataset" "$read_mode"
    expected=$(bigdata_expected_md5 "$dataset" "$read_mode") ||
        die "no big-data expected MD5 for $dataset $read_mode"
    rm -f "$BIGDATA_CASE_SAM" "$BIGDATA_CASE_LOG" "$BIGDATA_MD5_MARKER"
    build_bigdata_mem_command "$io_mode" "$dataset" "$read_mode" \
                              "$BIGDATA_CASE_SAM"

    echo
    echo "============================================================"
    echo "[BIG ] $dataset $read_mode"
    echo "[CONF] $configuration"
    echo "[I/O ] $io_mode"
    echo "[SAM ] $BIGDATA_CASE_SAM"
    echo "[LOG ] $BIGDATA_CASE_LOG"
    echo "============================================================"

    if ! submit_bigdata_job "$profile" "$job_name" "$BIGDATA_CASE_LOG"; then
        BIGDATA_VERIFY_STATUS=ERROR
        BIGDATA_VERIFY_EXPECTED=$expected
        BIGDATA_VERIFY_ACTUAL=$actual
        printf '[ERROR] %-18s %-2s  SWBWA job failed; partial SAM retained\n' \
               "$dataset" "$read_mode"
        return 1
    fi
    if [[ ! -f "$BIGDATA_CASE_SAM" ]]; then
        BIGDATA_VERIFY_STATUS=ERROR
        BIGDATA_VERIFY_EXPECTED=$expected
        BIGDATA_VERIFY_ACTUAL=$actual
        printf '[ERROR] %-18s %-2s  SAM output was not generated\n' \
               "$dataset" "$read_mode"
        return 1
    fi

    actual=$(file_md5 "$BIGDATA_CASE_SAM")
    BIGDATA_VERIFY_EXPECTED=$expected
    BIGDATA_VERIFY_ACTUAL=$actual
    if [[ "$actual" == "$expected" ]]; then
        rm -f "$BIGDATA_CASE_SAM"
        mkdir -p "$BIGDATA_STATUS_DIR"
        printf '%s\n' "$actual" > "$BIGDATA_MD5_MARKER"
        BIGDATA_VERIFY_STATUS=PASS
        printf '[PASS] %-18s %-2s  %s  (SAM removed)\n' \
               "$dataset" "$read_mode" "$actual"
        return 0
    fi

    BIGDATA_VERIFY_STATUS=FAIL
    printf '[FAIL] %-18s %-2s  expected=%s actual=%s  (SAM retained)\n' \
           "$dataset" "$read_mode" "$expected" "$actual"
    return 1
}

run_bigdata_configuration()
{
    local configuration=$1
    local profile=$2
    local exec_mode=$3
    local allocator=$4
    local config_dir="${BIGDATA_RESULT_ROOT}/${configuration}"
    local summary="${config_dir}/md5_results.tsv"
    local failures=0
    local io_mode
    local run_dir
    local dataset
    local read_mode

    echo
    echo "################################################################"
    echo "# BIG-DATA CONFIGURATION: $configuration"
    echo "################################################################"
    mkdir -p "$config_dir"
    if bigdata_configuration_requires_run "$config_dir"; then
        if ! build_configuration "$config_dir" "$exec_mode" "$allocator" \
                                 0 none ordered; then
            echo "[BUILD FAIL] $configuration" >&2
            return 1
        fi
    else
        echo "[BUILD SKIP] all selected big-data cases already passed"
    fi

    printf 'configuration\tio_mode\tdataset\tread_mode\tstatus\texpected\tactual\tlog\n' \
        > "$summary"
    for io_mode in "${IO_MODES[@]}"; do
        run_dir="${config_dir}/${io_mode}"
        mkdir -p "$run_dir"
        for dataset in "${SELECTED_DATASETS[@]}"; do
            for read_mode in PE SE; do
                bigdata_case_paths "$run_dir" "$dataset" "$read_mode"
                if bigdata_checkpoint_is_valid "$run_dir" "$dataset" \
                                                   "$read_mode"; then
                    BIGDATA_VERIFY_STATUS=PASS
                    BIGDATA_VERIFY_EXPECTED=$(
                        bigdata_expected_md5 "$dataset" "$read_mode")
                    read -r BIGDATA_VERIFY_ACTUAL < "$BIGDATA_MD5_MARKER"
                    printf '[SKIP] %-5s %-18s %-2s  MD5 already verified (%s)\n' \
                           "$io_mode" "$dataset" "$read_mode" \
                           "$BIGDATA_VERIFY_ACTUAL"
                elif ! run_bigdata_case "$configuration" "$profile" \
                                          "$io_mode" "$run_dir" "$dataset" \
                                          "$read_mode"; then
                    failures=$((failures + 1))
                fi
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                       "$configuration" "$io_mode" "$dataset" \
                       "$read_mode" "$BIGDATA_VERIFY_STATUS" \
                       "$BIGDATA_VERIFY_EXPECTED" \
                       "$BIGDATA_VERIFY_ACTUAL" "$BIGDATA_CASE_LOG" \
                       >> "$summary"
            done
        done
    done
    echo "[SUMMARY] $summary"
    ((failures == 0))
}

write_bigdata_global_summary()
{
    local output="${BIGDATA_RESULT_ROOT}/md5_results.tsv"
    local configuration
    local summary

    printf 'configuration\tio_mode\tdataset\tread_mode\tstatus\texpected\tactual\tlog\n' \
        > "$output"
    for configuration in single_system cgs_system cgs_pool \
                         cgs_cross_system cgs_cross_pool; do
        summary="${BIGDATA_RESULT_ROOT}/${configuration}/md5_results.tsv"
        [[ -f "$summary" ]] || continue
        tail -n +2 "$summary" >> "$output"
    done
    echo "[GLOBAL SUMMARY] $output"
}

bigdata_command()
{
    local datasets=all
    local failures=0

    while (($# > 0)); do
        case "$1" in
            --datasets) datasets=${2:?}; shift 2 ;;
            --result-root) BIGDATA_RESULT_ROOT=${2:?}; shift 2 ;;
            --chunk-bytes) CHUNK_BYTES=${2:?}; shift 2 ;;
            -h|--help) usage; exit 0 ;;
            *) die "unknown bigdata option '$1'" ;;
        esac
    done

    if [[ -n "$CHUNK_BYTES" ]]; then
        validate_positive_integer --chunk-bytes "$CHUNK_BYTES"
    fi
    validate_positive_integer BUILD_JOBS "$BUILD_JOBS"
    select_bigdata_datasets "$datasets"
    select_io_modes both
    [[ -d "$BIGDATA_DATA" ]] ||
        die "big-data input directory not found: $BIGDATA_DATA"
    [[ -f "$BIGDATA_MD5_FILE" ]] ||
        die "big-data MD5 reference not found: $BIGDATA_MD5_FILE"
    EXE=./SWBWA
    mkdir -p "$BIGDATA_RESULT_ROOT"

    if ! run_bigdata_configuration single_system single single system; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_configuration cgs_system cgs cgs system; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_configuration cgs_pool cgs cgs pool; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_configuration cgs_cross_system cgs_cross \
                                     cgs_cross system; then
        failures=$((failures + 1))
    fi
    if ! run_bigdata_configuration cgs_cross_pool cgs_cross \
                                     cgs_cross pool; then
        failures=$((failures + 1))
    fi
    write_bigdata_global_summary

    echo
    if ((failures == 0)); then
        echo "All big-data correctness configurations passed."
        return 0
    fi
    echo "$failures big-data configuration(s) had failures."
    return 1
}

main()
{
    local command=${1:-help}

    if (($# > 0)); then shift; fi
    case "$command" in
        run) run_command "$@" ;;
        verify) verify_command "$@" ;;
        matrix) matrix_command "$@" ;;
        bigdata) bigdata_command "$@" ;;
        bigdata-mpi) bigdata_mpi_command "$@" ;;
        help|-h|--help) usage ;;
        *) die "unknown command '$command'" ;;
    esac
}

main "$@"
