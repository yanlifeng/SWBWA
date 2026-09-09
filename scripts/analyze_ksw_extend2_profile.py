#!/usr/bin/env python3
"""Summarize a ksw_extend2 profile and estimate 8-way packing potential."""

import argparse
import statistics


def read_profile(path):
    meta = {}
    samples = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("# "):
                fields = line[2:].split("\t", 1)
                if len(fields) == 2:
                    meta[fields[0]] = fields[1]
                continue
            if line.startswith("sample_index\t"):
                continue
            fields = line.split("\t")
            if len(fields) != 8:
                continue
            samples.append({
                "qlen": int(fields[1]),
                "tlen": int(fields[2]),
                "requested_w": int(fields[3]),
                "effective_w": int(fields[4]),
                "rows": int(fields[5]),
                "dp_cells": int(fields[6]),
                "max_score": int(fields[7]),
            })
    return meta, samples


def quantile(values, fraction):
    if not values:
        return 0
    values = sorted(values)
    index = min(len(values) - 1, int(fraction * len(values)))
    return values[index]


def pack_stats(samples, width=8):
    groups = [samples[i:i + width] for i in range(0, len(samples), width)]
    groups = [group for group in groups if len(group) == width]
    if not groups:
        return None
    ratios = []
    utilization = []
    for group in groups:
        work = [item["dp_cells"] for item in group]
        largest = max(work)
        smallest = min(work)
        ratios.append(float(largest) / max(1, smallest))
        utilization.append(float(sum(work)) / (width * largest))
    return {
        "groups": len(groups),
        "ratio_median": statistics.median(ratios),
        "ratio_p90": quantile(ratios, 0.90),
        "util_median": statistics.median(utilization),
        "util_p10": quantile(utilization, 0.10),
        "ratio_le_1_25": sum(r <= 1.25 for r in ratios) / float(len(ratios)),
        "ratio_le_1_50": sum(r <= 1.50 for r in ratios) / float(len(ratios)),
    }


def print_summary(path, meta, samples):
    calls = int(meta.get("calls", 0))
    total_cells = int(meta.get("total_dp_cells", 0))
    total_rows = int(meta.get("total_rows", 0))
    print("# ksw_extend2 distribution: %s" % path)
    print("calls=%d samples=%d sample_fraction=%.6f" % (
        calls, len(samples), len(samples) / float(max(1, calls))))
    print("total_dp_cells=%d total_rows=%d avg_cells_per_call=%.2f" % (
        total_cells, total_rows, total_cells / float(max(1, calls))))
    if not samples:
        return
    for field in ("qlen", "tlen", "requested_w", "effective_w", "rows", "dp_cells"):
        values = [item[field] for item in samples]
        print("%-13s min=%-8d p50=%-8d p90=%-8d p99=%-8d max=%-8d mean=%.2f" % (
            field, min(values), quantile(values, 0.50), quantile(values, 0.90),
            quantile(values, 0.99), max(values), statistics.mean(values)))

    stats = pack_stats(samples)
    if stats is not None:
        print("pack8 consecutive sampled calls: groups=%d max/min median=%.3f p90=%.3f" % (
            stats["groups"], stats["ratio_median"], stats["ratio_p90"]))
        print("pack8 utilization: median=%.3f p10=%.3f; groups max/min<=1.25: %.1f%%; <=1.50: %.1f%%" % (
            stats["util_median"], stats["util_p10"],
            100.0 * stats["ratio_le_1_25"], 100.0 * stats["ratio_le_1_50"]))

        sorted_samples = sorted(samples, key=lambda item: item["dp_cells"])
        sorted_stats = pack_stats(sorted_samples)
        print("pack8 after ideal offline work sorting: utilization median=%.3f max/min median=%.3f" % (
            sorted_stats["util_median"], sorted_stats["ratio_median"]))

    print("\nInterpretation:")
    print("- dp_cells is the actual number of band cells visited by the scalar extension call.")
    print("- consecutive pack8 is a no-reordering baseline; samples are recorded in multi-thread completion order.")
    print("- sorted pack8 is an upper bound; it is not a valid online scheduling strategy by itself.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", nargs="+")
    args = parser.parse_args()
    for path in args.profile:
        meta, samples = read_profile(path)
        print_summary(path, meta, samples)
        if path != args.profile[-1]:
            print("\n" + "=" * 72 + "\n")


if __name__ == "__main__":
    main()
