#!/usr/bin/env bash

set -euo pipefail

usage()
{
    cat >&2 <<'EOF'
Usage:
  ./run.sh single [--] SWBWA mem arguments...
  ./run.sh cgs [--] SWBWA mem arguments...
  ./run.sh cgs_cross [--] SWBWA mem arguments...
  ./run.sh mpi [--nodes N] [--ranks N] [--] SWBWA mem arguments...

The executable must already have been built with the matching Makefile mode.
This wrapper only assembles the Sunway bsub request; it does not rebuild.
EOF
    exit 2
}

[[ $# -gt 0 ]] || usage
profile=$1
shift
nodes=1
ranks=1
queue=${QUEUE:-q_share}

while [[ $# -gt 0 ]]; do
    case $1 in
        --nodes) nodes=$2; shift 2 ;;
        --ranks) ranks=$2; shift 2 ;;
        --queue) queue=$2; shift 2 ;;
        --) shift; break ;;
        *) break ;;
    esac
done
[[ $# -gt 0 ]] || usage

case $profile in
    single)
        bsub_args=(-n 1 -cgsp 64 -share_size "${SHARE_SIZE:-12000}") ;;
    cgs|cgs_cross)
        bsub_args=(-n 1 -cgsp 64 -mpecg 6
                   -share_size "${SHARE_SIZE:-2000}"
                   -xmalloc -cross_size "${CROSS_SIZE:-42000}") ;;
    mpi)
        bsub_args=(-N "$nodes" -np "$ranks" -cgsp 64
                   -share_size "${SHARE_SIZE:-12000}") ;;
    *) usage ;;
esac

exec bsub -I -b -q "$queue" "${bsub_args[@]}" "$@"
