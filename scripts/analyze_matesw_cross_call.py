#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


PATTERNS = {
    "candidates": re.compile(r"valid KSW candidates:\s+(\d+)"),
    "sam_calls": re.compile(r"mem_sam_pe calls:\s+(\d+)"),
    "both_directions": re.compile(
        r"both mate directions active:\s+(\d+) \(([0-9.]+)%\)"
    ),
    "same_pe": re.compile(
        r"same-PE direction pairs:\s+(\d+) "
        r"\((\d+) candidates, ([0-9.]+)% coverage\)"
    ),
    "adjacent_pe": re.compile(
        r"adjacent-PE pairs:\s+(\d+) "
        r"\((\d+) candidates, ([0-9.]+)% coverage\)"
    ),
}


def match(pattern, text, label):
    result = pattern.search(text)
    if result is None:
        raise ValueError(f"missing {label}")
    return result.groups()


def ideal_speedup(candidates, pairs):
    calls = candidates - pairs
    return candidates / calls if calls else 1.0


def parse_log(path):
    text = path.read_text(errors="replace")
    candidates = int(match(PATTERNS["candidates"], text, "candidates")[0])
    sam_calls = int(match(PATTERNS["sam_calls"], text, "sam_calls")[0])
    both_count, both_pct = match(
        PATTERNS["both_directions"], text, "both directions"
    )
    same_pairs, same_candidates, same_pct = match(
        PATTERNS["same_pe"], text, "same-PE pairs"
    )
    adjacent_pairs, adjacent_candidates, adjacent_pct = match(
        PATTERNS["adjacent_pe"], text, "adjacent-PE pairs"
    )
    return {
        "dataset": path.stem,
        "sam_calls": sam_calls,
        "candidates": candidates,
        "both_count": int(both_count),
        "both_pct": float(both_pct),
        "same_pairs": int(same_pairs),
        "same_candidates": int(same_candidates),
        "same_pct": float(same_pct),
        "same_speedup": ideal_speedup(candidates, int(same_pairs)),
        "adjacent_pairs": int(adjacent_pairs),
        "adjacent_candidates": int(adjacent_candidates),
        "adjacent_pct": float(adjacent_pct),
        "adjacent_speedup": ideal_speedup(candidates, int(adjacent_pairs)),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Summarize cross-call mate-SW pairing opportunities."
    )
    parser.add_argument("logs", nargs="+", type=Path)
    args = parser.parse_args()

    print(
        "dataset\tsam_pe_calls\tksw_candidates\tboth_directions_pct\t"
        "same_pe_coverage_pct\tsame_pe_ideal_speedup\t"
        "adjacent_pe_coverage_pct\tadjacent_pe_ideal_speedup"
    )
    for path in args.logs:
        row = parse_log(path)
        print(
            f"{row['dataset']}\t{row['sam_calls']}\t{row['candidates']}\t"
            f"{row['both_pct']:.2f}\t{row['same_pct']:.2f}\t"
            f"{row['same_speedup']:.3f}\t{row['adjacent_pct']:.2f}\t"
            f"{row['adjacent_speedup']:.3f}"
        )


if __name__ == "__main__":
    main()
