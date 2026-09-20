#!/usr/bin/env bash
# Direct SSH, logs only, no deletion; run after the experiment has stopped.
set -euo pipefail

REMOTE=${REMOTE:?Set REMOTE to the SSH destination}
SSH_PORT=${SSH_PORT:-22}
AB_REMOTE=${AB_REMOTE:?Set AB_REMOTE to the remote experiment directory}
LOCAL=${LOCAL:-$(cd "$(dirname "$0")" && pwd)/fetch_raw}
[[ "$AB_REMOTE" =~ ^/[a-zA-Z0-9_./-]+$ ]] || exit 2

filters=(--include='*/' --include='*.log' --include='*.out'
         --include='*.joblog' --include='*.status' --include='*.tsv' --exclude='*')
mkdir -p "$LOCAL"
stats=$(LC_ALL=C rsync -an --stats "${filters[@]}" -e "ssh -p $SSH_PORT" \
    "$REMOTE:$AB_REMOTE/logs/" "$LOCAL/")
bytes=$(printf '%s\n' "$stats" | awk '/^Total file size:/ {gsub(/,/, "", $4); print $4}')
[[ "$bytes" =~ ^[0-9]+$ && "$bytes" -le 10485760 ]] || {
    echo "Refusing transfer: selected logs exceed 10 MiB or size is unknown" >&2
    exit 1
}
rsync -av --bwlimit=500 "${filters[@]}" -e "ssh -p $SSH_PORT" \
    "$REMOTE:$AB_REMOTE/logs/" "$LOCAL/"
echo "fetched logs to $LOCAL (preflight: $bytes bytes)"
