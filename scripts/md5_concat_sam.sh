#!/usr/bin/env bash
set -euo pipefail

if (($# < 2)); then
    echo "Usage: $0 RESULT INPUT..." >&2
    exit 2
fi

result=$1
shift

if command -v md5sum >/dev/null 2>&1; then
    cat "$@" | md5sum | awk '{print $1}' > "$result"
elif command -v md5 >/dev/null 2>&1; then
    cat "$@" | md5 > "$result"
else
    echo "Neither md5sum nor md5 is available" >&2
    exit 1
fi
