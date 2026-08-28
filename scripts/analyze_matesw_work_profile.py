#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def required(pattern, text, label):
    match = re.search(pattern, text, re.MULTILINE)
    if match is None:
        raise ValueError(f"missing {label}")
    return match


def scaled_number(value):
    suffixes = {"K": 1e3, "M": 1e6, "G": 1e9}
    suffix = value[-1]
    if suffix in suffixes:
        return float(value[:-1]) * suffixes[suffix]
    return float(value)


def kernel_cycles(text, kernel):
    match = required(
        rf"^\|\s*{re.escape(kernel)}\s*\|\s*([0-9.]+[KMG]?)\|",
        text,
        f"{kernel} cycles",
    )
    return scaled_number(match.group(1))


def parse_log(path):
    text = path.read_text(errors="replace")

    profiled = int(required(
        r"profiled ordered pairs:\s+(\d+)", text, "profiled pairs"
    ).group(1))

    def percent(label):
        return float(required(
            rf"{re.escape(label)}:\s+\d+ \(([0-9.]+)%\)",
            text,
            label,
        ).group(1))

    def ratio_percent(label):
        return float(required(
            rf"{re.escape(label)}:\s+<=1\.10=\d+ \(([0-9.]+)%\)",
            text,
            label,
        ).group(1))

    forward = required(
        r"all forward:\s+(\d+) \+ (\d+) \(Lazy-F ([0-9.]+)%\)",
        text,
        "forward work",
    )
    reverse = required(
        r"all reverse:\s+(\d+) \+ (\d+) \(Lazy-F ([0-9.]+)%\)",
        text,
        "reverse work",
    )
    paired_forward = required(
        r"paired forward:\s+serial=\d+ lockstep=\d+ \(([0-9.]+)x\)",
        text,
        "paired forward work",
    )
    paired_reverse = required(
        r"paired reverse:\s+serial=\d+ lockstep=\d+ \(([0-9.]+)x\)",
        text,
        "paired reverse work",
    )
    projected = required(
        r"all candidate DP projected:\s+original=\d+ lockstep=\d+ "
        r"\(([0-9.]+)x\)",
        text,
        "projected DP work",
    )
    stage2 = float(required(
        r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s",
        text,
        "stage 2 time",
    ).group(1))
    part3 = float(required(
        r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s",
        text,
        "part 3 time",
    ).group(1))

    worker_cycles = kernel_cycles(text, "WORKER ALIGNMENT")
    dp_cycles = (kernel_cycles(text, "KSW DP FORWARD") +
                 kernel_cycles(text, "KSW DP REVERSE"))
    dp_speedup = float(projected.group(1))
    dp_cycle_share = dp_cycles / worker_cycles
    stage2_upper = 1.0 / (
        1.0 - dp_cycle_share * (1.0 - 1.0 / dp_speedup)
    )

    return {
        "dataset": path.stem,
        "stage2": stage2,
        "part3": part3,
        "profiled_pairs": profiled,
        "exact_qlen_pct": percent("exact forward qlen pairs"),
        "exact_tlen_pct": percent("exact forward tlen pairs"),
        "exact_forward_dimensions_pct": percent("exact forward dimensions"),
        "exact_reverse_dimensions_pct": percent("exact reverse dimensions"),
        "exact_dimensions_pct": percent("exact both-pass dimensions"),
        "qlen_110_pct": ratio_percent("qlen ratio max/min"),
        "tlen_110_pct": ratio_percent("tlen ratio max/min"),
        "forward_110_pct": ratio_percent("forward work ratio max/min"),
        "reverse_110_pct": ratio_percent("reverse work ratio max/min"),
        "total_110_pct": ratio_percent("total DP work ratio max/min"),
        "forward_lazy_pct": float(forward.group(3)),
        "reverse_lazy_pct": float(reverse.group(3)),
        "paired_forward_speedup": float(paired_forward.group(1)),
        "paired_reverse_speedup": float(paired_reverse.group(1)),
        "projected_dp_speedup": dp_speedup,
        "dp_cycle_share_pct": 100.0 * dp_cycle_share,
        "stage2_upper": stage2_upper,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Summarize same-PE mate-SW DP-work compatibility."
    )
    parser.add_argument("logs", nargs="+", type=Path)
    args = parser.parse_args()

    columns = [
        "dataset", "stage2_s", "part3_s", "profiled_pairs",
        "exact_qlen_pct", "exact_tlen_pct", "exact_forward_dimensions_pct",
        "exact_reverse_dimensions_pct", "exact_both_pass_dimensions_pct",
        "qlen_le_1.10_pct", "tlen_le_1.10_pct",
        "forward_work_le_1.10_pct", "reverse_work_le_1.10_pct",
        "total_work_le_1.10_pct", "forward_lazy_pct",
        "reverse_lazy_pct", "paired_forward_speedup",
        "paired_reverse_speedup", "projected_all_dp_speedup",
        "dp_cycle_share_pct", "projected_stage2_upper",
    ]
    print("\t".join(columns))
    for path in args.logs:
        row = parse_log(path)
        values = [
            row["dataset"], f"{row['stage2']:.3f}", f"{row['part3']:.3f}",
            str(row["profiled_pairs"]), f"{row['exact_qlen_pct']:.2f}",
            f"{row['exact_tlen_pct']:.2f}",
            f"{row['exact_forward_dimensions_pct']:.2f}",
            f"{row['exact_reverse_dimensions_pct']:.2f}",
            f"{row['exact_dimensions_pct']:.2f}",
            f"{row['qlen_110_pct']:.2f}", f"{row['tlen_110_pct']:.2f}",
            f"{row['forward_110_pct']:.2f}",
            f"{row['reverse_110_pct']:.2f}", f"{row['total_110_pct']:.2f}",
            f"{row['forward_lazy_pct']:.2f}",
            f"{row['reverse_lazy_pct']:.2f}",
            f"{row['paired_forward_speedup']:.3f}",
            f"{row['paired_reverse_speedup']:.3f}",
            f"{row['projected_dp_speedup']:.3f}",
            f"{row['dp_cycle_share_pct']:.2f}", f"{row['stage2_upper']:.3f}",
        ]
        print("\t".join(values))


if __name__ == "__main__":
    main()
