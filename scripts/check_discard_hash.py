#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


HASH_RE = re.compile(
    r"\[SWBWA output hash rank\s+(\d+)/(\d+)\] "
    r"calls=(\d+) bytes=(\d+) sum=0x([0-9a-fA-F]{16}) "
    r"xor=0x([0-9a-fA-F]{16}) enabled=([01])"
)
MASK64 = (1 << 64) - 1


def parse_log(path):
    ranks = {}
    world_size = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = HASH_RE.search(line)
            if match is None:
                continue
            rank = int(match.group(1))
            size = int(match.group(2))
            enabled = int(match.group(7))
            if size <= 0 or rank < 0 or rank >= size:
                raise ValueError(
                    f"{path}: invalid rank/world size {rank}/{size}"
                )
            if world_size is not None and world_size != size:
                raise ValueError(f"{path}: inconsistent MPI world sizes")
            if rank in ranks:
                raise ValueError(f"{path}: duplicate hash for rank {rank}")
            if not enabled:
                raise ValueError(f"{path}: hash disabled for rank {rank}")
            world_size = size
            ranks[rank] = {
                "calls": int(match.group(3)),
                "bytes": int(match.group(4)),
                "sum": int(match.group(5), 16),
                "xor": int(match.group(6), 16),
            }

    if world_size is None:
        raise ValueError(f"{path}: no SWBWA discard hash lines found")
    missing = sorted(set(range(world_size)) - set(ranks))
    if missing:
        raise ValueError(f"{path}: missing ranks {missing}")

    calls = sum(item["calls"] for item in ranks.values())
    byte_count = sum(item["bytes"] for item in ranks.values())
    hash_sum = sum(item["sum"] for item in ranks.values()) & MASK64
    hash_xor = 0
    for item in ranks.values():
        hash_xor ^= item["xor"]
    return world_size, calls, byte_count, hash_sum, hash_xor


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-rank SWBWA discard hashes without MPI."
    )
    parser.add_argument(
        "--require-matching",
        action="store_true",
        help="fail unless all input logs have the same aggregate fingerprint",
    )
    parser.add_argument("logs", nargs="+", type=Path)
    args = parser.parse_args()

    print("log\tranks\tcalls\tbytes\thash_sum\thash_xor")
    failed = False
    fingerprints = []
    for path in args.logs:
        try:
            size, calls, byte_count, hash_sum, hash_xor = parse_log(path)
        except (OSError, ValueError) as error:
            print(f"[ERROR] {error}", file=sys.stderr)
            failed = True
            continue
        print(
            f"{path}\t{size}\t{calls}\t{byte_count}\t"
            f"{hash_sum:016x}\t{hash_xor:016x}"
        )
        fingerprints.append((calls, byte_count, hash_sum, hash_xor))
    if args.require_matching and fingerprints and any(
        item != fingerprints[0] for item in fingerprints[1:]
    ):
        print("[ERROR] aggregate fingerprints do not match", file=sys.stderr)
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
