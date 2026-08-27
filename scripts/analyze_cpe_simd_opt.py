#!/usr/bin/env python3
"""Summarize the CPE SIMD/BWT experiment and render dependency-free SVGs."""

import argparse
import html
import re
from pathlib import Path


PROFILE_REGIONS = [
    "WORKER_ALIGNMENT",
    "MEM_CHAIN",
    "CHAIN_FILTER",
    "CHAIN_EXTENSION",
    "ALIGNMENT_FINALIZE",
    "MATE_RESCUE",
    "PAIRING",
    "SAM_FORMAT",
    "SAM_COPY",
    "MEM_CHAIN_COLLECT",
    "MEM_CHAIN_BUILD",
    "CHAIN_EXTENSION_DP",
    "MATE_REF_FETCH",
    "MATE_KSW_ALIGN",
    "KSW_QUERY_INIT_FORWARD",
    "KSW_DP_FORWARD",
    "KSW_QUERY_INIT_REVERSE",
    "KSW_DP_REVERSE",
    "MATE_DEDUP",
    "DEDUP_SORT_END",
    "DEDUP_REDUNDANCY",
    "DEDUP_SORT_SCORE",
    "MEM_COLLECT_FIRST",
    "MEM_COLLECT_SPLIT",
    "MEM_COLLECT_LAST",
    "MEM_COLLECT_SORT",
    "CHAIN_BUILD_REPETITIVE",
    "CHAIN_BUILD_SA",
    "CHAIN_BUILD_RID",
    "CHAIN_BUILD_TREE_SEARCH",
    "CHAIN_BUILD_MERGE",
    "CHAIN_BUILD_INSERT",
    "CHAIN_BUILD_FINALIZE",
]

PARENT = {
    "WORKER_ALIGNMENT": None,
    "MEM_CHAIN": "WORKER_ALIGNMENT",
    "MEM_CHAIN_COLLECT": "MEM_CHAIN",
    "MEM_COLLECT_FIRST": "MEM_CHAIN_COLLECT",
    "MEM_COLLECT_SPLIT": "MEM_CHAIN_COLLECT",
    "MEM_COLLECT_LAST": "MEM_CHAIN_COLLECT",
    "MEM_COLLECT_SORT": "MEM_CHAIN_COLLECT",
    "MEM_CHAIN_BUILD": "MEM_CHAIN",
    "CHAIN_BUILD_REPETITIVE": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_SA": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_RID": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_TREE_SEARCH": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_MERGE": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_INSERT": "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_FINALIZE": "MEM_CHAIN_BUILD",
    "CHAIN_FILTER": "WORKER_ALIGNMENT",
    "CHAIN_EXTENSION": "WORKER_ALIGNMENT",
    "CHAIN_EXTENSION_DP": "CHAIN_EXTENSION",
    "ALIGNMENT_FINALIZE": "WORKER_ALIGNMENT",
    "SAM_FORMAT": "WORKER_ALIGNMENT",
    "MATE_RESCUE": "SAM_FORMAT",
    "MATE_REF_FETCH": "MATE_RESCUE",
    "MATE_KSW_ALIGN": "MATE_RESCUE",
    "KSW_QUERY_INIT_FORWARD": "MATE_KSW_ALIGN",
    "KSW_DP_FORWARD": "MATE_KSW_ALIGN",
    "KSW_QUERY_INIT_REVERSE": "MATE_KSW_ALIGN",
    "KSW_DP_REVERSE": "MATE_KSW_ALIGN",
    "MATE_DEDUP": "MATE_RESCUE",
    "DEDUP_SORT_END": "MATE_DEDUP",
    "DEDUP_REDUNDANCY": "MATE_DEDUP",
    "DEDUP_SORT_SCORE": "MATE_DEDUP",
    "PAIRING": "SAM_FORMAT",
    "SAM_COPY": None,
}

TREE_ORDER = [
    "WORKER_ALIGNMENT",
    "MEM_CHAIN",
    "MEM_CHAIN_COLLECT",
    "MEM_COLLECT_FIRST",
    "MEM_COLLECT_SPLIT",
    "MEM_COLLECT_LAST",
    "MEM_COLLECT_SORT",
    "MEM_CHAIN_BUILD",
    "CHAIN_BUILD_REPETITIVE",
    "CHAIN_BUILD_SA",
    "CHAIN_BUILD_RID",
    "CHAIN_BUILD_TREE_SEARCH",
    "CHAIN_BUILD_MERGE",
    "CHAIN_BUILD_INSERT",
    "CHAIN_BUILD_FINALIZE",
    "CHAIN_FILTER",
    "CHAIN_EXTENSION",
    "CHAIN_EXTENSION_DP",
    "ALIGNMENT_FINALIZE",
    "SAM_FORMAT",
    "MATE_RESCUE",
    "MATE_REF_FETCH",
    "MATE_KSW_ALIGN",
    "KSW_QUERY_INIT_FORWARD",
    "KSW_DP_FORWARD",
    "KSW_QUERY_INIT_REVERSE",
    "KSW_DP_REVERSE",
    "MATE_DEDUP",
    "DEDUP_SORT_END",
    "DEDUP_REDUNDANCY",
    "DEDUP_SORT_SCORE",
    "PAIRING",
]

RUNS = [
    ("int32", "ERR 75 bp", "int32 baseline", "int32_baseline/ERR1203383_PE.run.log", True),
    ("fp16x32", "ERR 75 bp", "FP16 32 lanes (rejected)", "fp16/ERR1203383_PE.run.log", False),
    ("fp16x16", "ERR 75 bp", "FP16 16 logical lanes", "fp16x16/ERR1203383_PE.run.log", True),
    ("sa4", "ERR 75 bp", "SA batch 4", "sa_batch4/ERR1203383_PE.run.log", True),
    ("sa8", "ERR 75 bp", "SA batch 8", "sa_batch8/ERR1203383_PE.run.log", True),
    ("sa16", "ERR 75 bp", "SA batch 16", "sa_batch16/ERR1203383_PE.run.log", True),
    ("rid_ldm", "ERR 75 bp", "contig offsets in LDM", "rid_ldm/ERR1203383_PE.run.log", True),
    ("combined", "ERR 75 bp", "combined before allocation fix", "chain_arena/ERR1203383_PE.run.log", True),
    ("final", "ERR 75 bp", "final compact KSW allocation", "final_profile1/ERR1203383_PE.run.log", True),
    ("prefetch8", "ERR 75 bp", "lookahead 8 (rejected)", "prefetch_distance8/ERR1203383_PE.run.log", True),
    ("fp16x32", "SRR 150 bp", "FP16 32 lanes (rejected)", "fp16/SRR7963242_PE.run.log", False),
    ("fp16x16", "SRR 150 bp", "FP16 16 logical lanes", "fp16x16/SRR7963242_PE.run.log", True),
    ("combined", "SRR 150 bp", "combined before allocation fix", "chain_arena/SRR7963242_PE.run.log", True),
    ("final", "SRR 150 bp", "final compact KSW allocation", "final_profile1/SRR7963242_PE.run.log", True),
]

NUMBER_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]*)?)([KMGT]?)\s*$")


def scaled_number(value):
    match = NUMBER_RE.match(value)
    if match is None:
        raise ValueError("invalid LWPF value: {!r}".format(value))
    scale = {"": 1.0, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    return float(match.group(1)) * scale[match.group(2)]


def last_float(pattern, text):
    values = re.findall(pattern, text, flags=re.MULTILINE)
    if not values:
        raise ValueError("missing timing pattern: {}".format(pattern))
    return float(values[-1])


def parse_timing(path):
    text = path.read_text(errors="replace")
    return {
        "stage2": last_float(
            r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s", text),
        "part3": last_float(
            r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s", text),
    }


def parse_lwpf(path):
    text = path.read_text(errors="replace")
    rows = []
    in_table = False
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\r")
        if line.startswith("LWPF kernel summary"):
            rows = []
            in_table = True
            continue
        if in_table and line.startswith("===="):
            break
        if not in_table or not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.split("|")[1:-1]]
        if len(cells) != 25:
            continue
        rows.append([scaled_number(cell) for cell in cells[1:]])
    if len(rows) not in (22, 33):
        raise ValueError("{}: unexpected LWPF row count {}".format(path, len(rows)))
    return {
        name: {
            "cycles_avg": values[0],
            "cycles_min": values[1],
            "cycles_max": values[2],
            "instructions_avg": values[3],
            "dcache_accesses_avg": values[12],
            "dcache_misses_avg": values[15],
            "instruction_stall_cycles_avg": values[21],
        }
        for name, values in zip(PROFILE_REGIONS, rows)
    }


def parse_dedup_hist(path):
    text = path.read_text(errors="replace")
    match = re.search(
        r"Mate dedup input-size distribution \(([0-9]+) calls\)(.*?)(?:={20,})",
        text, flags=re.DOTALL)
    if match is None:
        raise ValueError("{}: missing dedup histogram".format(path))
    rows = []
    for label, count, percent in re.findall(
            r"^\s*([0-9]+(?:-[0-9]+)?|65\+)\s*:\s*([0-9]+)\s+([0-9.]+)%",
            match.group(2), flags=re.MULTILINE):
        rows.append((label, int(count), float(percent)))
    return int(match.group(1)), rows


def depth(region):
    value = 0
    parent = PARENT[region]
    while parent is not None:
        value += 1
        parent = PARENT[parent]
    return value


def svg_header(width, height, title, subtitle):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
        'viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text { font-family: Arial, sans-serif; fill: #18272f; '
        'letter-spacing: 0; } .small { font-size: 12px; fill: #5d6870; }</style>',
        '<text x="32" y="38" font-size="21" font-weight="700">{}</text>'.format(
            html.escape(title)),
        '<text x="32" y="60" class="small">{}</text>'.format(
            html.escape(subtitle)),
    ]


def write_run_summary(rows, output):
    lines = ["dataset\tvariant\tlabel\tcorrect\tstage2_s\tpart3_s"]
    for row in rows:
        lines.append("\t".join([
            row["dataset"], row["variant"], row["label"],
            "yes" if row["correct"] else "no", "{:.3f}".format(row["stage2"]),
            "{:.3f}".format(row["part3"]),
        ]))
    output.write_text("\n".join(lines) + "\n")


def write_stage2_svg(rows, output):
    width = 1280
    row_height = 33
    height = 105 + row_height * len(rows)
    left = 330
    maximum = max(row["stage2"] for row in rows)
    out = svg_header(
        width, height, "CPE SIMD and BWT optimization experiments",
        "Stage 2 wall time; hatched red variants changed output and were rejected",
    )
    previous_dataset = None
    for index, row in enumerate(rows):
        y = 82 + index * row_height
        if row["dataset"] != previous_dataset:
            out.append('<line x1="24" y1="{}" x2="1250" y2="{}" stroke="#d7dcdf"/>'.format(
                y - 16, y - 16))
            previous_dataset = row["dataset"]
        label = "{} / {}".format(row["dataset"], row["label"])
        length = 760 * row["stage2"] / maximum
        color = "#087e8b" if row["correct"] else "#c44e52"
        out.append('<text x="32" y="{}" font-size="12">{}</text>'.format(
            y + 8, html.escape(label)))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="12" fill="{}"/>'.format(
            left, y - 3, length, color))
        out.append('<text x="{}" y="{}" class="small">{:.3f} s; part3 {:.3f} s</text>'.format(
            left + length + 8, y + 8, row["stage2"], row["part3"]))
    out.append("</svg>")
    output.write_text("\n".join(out) + "\n")


def write_hotspots(profiles, tsv_path, svg_path):
    lines = [
        "dataset\tregion\tparent\tdepth\tcycles_avg\tcycles_min\tcycles_max\t"
        "parent_percent\tworker_percent\tinstructions_avg\tdcache_misses_avg\t"
        "instruction_stall_cycles_avg"
    ]
    for dataset, profile in profiles.items():
        worker_cycles = profile["WORKER_ALIGNMENT"]["cycles_avg"]
        for region in TREE_ORDER:
            if region not in profile:
                continue
            parent = PARENT[region]
            parent_cycles = profile[parent]["cycles_avg"] if parent else worker_cycles
            values = profile[region]
            lines.append("\t".join([
                dataset, region, parent or "", str(depth(region)),
                "{:.0f}".format(values["cycles_avg"]),
                "{:.0f}".format(values["cycles_min"]),
                "{:.0f}".format(values["cycles_max"]),
                "{:.3f}".format(100.0 * values["cycles_avg"] / parent_cycles),
                "{:.3f}".format(100.0 * values["cycles_avg"] / worker_cycles),
                "{:.0f}".format(values["instructions_avg"]),
                "{:.0f}".format(values["dcache_misses_avg"]),
                "{:.0f}".format(values["instruction_stall_cycles_avg"]),
            ]))
    tsv_path.write_text("\n".join(lines) + "\n")

    width = 1500
    row_height = 29
    height = 105 + row_height * len(TREE_ORDER)
    bar_x = 390
    colors = {"ERR 75 bp": "#0072b2", "SRR 150 bp": "#d55e00"}
    out = svg_header(
        width, height, "Final CPE hotspot hierarchy",
        "Percent of WORKER_ALIGNMENT cycles; indented rows are children, not additive siblings",
    )
    for legend_index, dataset in enumerate(profiles):
        x = 1040 + legend_index * 160
        out.append('<rect x="{}" y="29" width="13" height="9" fill="{}"/>'.format(
            x, colors[dataset]))
        out.append('<text x="{}" y="38" class="small">{}</text>'.format(
            x + 18, html.escape(dataset)))
    for index, region in enumerate(TREE_ORDER):
        y = 79 + index * row_height
        indent = depth(region)
        out.append('<text x="{}" y="{}" font-size="12">{}</text>'.format(
            32 + indent * 14, y + 8, html.escape(region.replace("_", " ").title())))
        for profile_index, (dataset, profile) in enumerate(profiles.items()):
            if region not in profile:
                continue
            percent = 100.0 * profile[region]["cycles_avg"] / profile["WORKER_ALIGNMENT"]["cycles_avg"]
            length = 820 * percent / 100.0
            bar_y = y - 4 + profile_index * 9
            out.append('<rect x="{}" y="{}" width="{:.2f}" height="7" fill="{}"/>'.format(
                bar_x, bar_y, length, colors[dataset]))
            out.append('<text x="{}" y="{}" class="small">{:.1f}%</text>'.format(
                bar_x + length + 6, bar_y + 7, percent))
    out.append("</svg>")
    svg_path.write_text("\n".join(out) + "\n")


def write_dedup(root, tsv_path, svg_path):
    datasets = {
        "ERR 75 bp": root / "dedup_hist/ERR1203383_PE.run.log",
        "SRR 150 bp": root / "dedup_hist/SRR7963242_PE.run.log",
    }
    parsed = {name: parse_dedup_hist(path) for name, path in datasets.items()}
    lines = ["dataset\ttotal_calls\tbucket\tcount\tpercent"]
    for dataset, (total, rows) in parsed.items():
        for label, count, percent in rows:
            lines.append("\t".join([
                dataset, str(total), label, str(count), "{:.2f}".format(percent)]))
    tsv_path.write_text("\n".join(lines) + "\n")

    labels = [row[0] for row in next(iter(parsed.values()))[1]]
    width, height = 1250, 500
    out = svg_header(
        width, height, "Mate dedup input-size distribution",
        "Most calls exceed 16 candidates, so a tiny-array-only sorter has limited coverage",
    )
    chart_left, chart_top, chart_height = 75, 90, 320
    group_width = 1040 / len(labels)
    colors = {"ERR 75 bp": "#0072b2", "SRR 150 bp": "#d55e00"}
    for tick in range(0, 51, 10):
        y = chart_top + chart_height * (1.0 - tick / 55.0)
        out.append('<line x1="{}" y1="{}" x2="1170" y2="{}" stroke="#dde1e3"/>'.format(
            chart_left, y, y))
        out.append('<text x="30" y="{}" class="small">{}%</text>'.format(y + 4, tick))
    for index, label in enumerate(labels):
        x = chart_left + index * group_width
        out.append('<text x="{}" y="435" class="small">{}</text>'.format(
            x + 8, html.escape(label)))
        for dataset_index, (dataset, (_, rows)) in enumerate(parsed.items()):
            percent = rows[index][2]
            bar_height = chart_height * percent / 55.0
            bar_x = x + 8 + dataset_index * 25
            out.append('<rect x="{}" y="{:.2f}" width="20" height="{:.2f}" fill="{}"/>'.format(
                bar_x, chart_top + chart_height - bar_height, bar_height, colors[dataset]))
    for index, dataset in enumerate(parsed):
        x = 880 + index * 150
        out.append('<rect x="{}" y="29" width="13" height="9" fill="{}"/>'.format(
            x, colors[dataset]))
        out.append('<text x="{}" y="38" class="small">{}</text>'.format(
            x + 18, html.escape(dataset)))
    out.append("</svg>")
    svg_path.write_text("\n".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path,
        default=Path("correctness_results/cpe_simd_opt_20260826"))
    args = parser.parse_args()
    output = args.root / "analysis"
    output.mkdir(parents=True, exist_ok=True)

    rows = []
    for variant, dataset, label, relative, correct in RUNS:
        timing = parse_timing(args.root / relative)
        rows.append({
            "variant": variant, "dataset": dataset, "label": label,
            "correct": correct, **timing,
        })
    write_run_summary(rows, output / "run_summary.tsv")
    write_stage2_svg(rows, output / "stage2_comparison.svg")

    profiles = {
        "ERR 75 bp": parse_lwpf(args.root / "final_profile1/ERR1203383_PE.run.log"),
        "SRR 150 bp": parse_lwpf(args.root / "final_profile1/SRR7963242_PE.run.log"),
    }
    write_hotspots(
        profiles, output / "hotspot_hierarchy.tsv",
        output / "hotspot_hierarchy.svg")
    write_dedup(
        args.root, output / "mate_dedup_distribution.tsv",
        output / "mate_dedup_distribution.svg")


if __name__ == "__main__":
    main()
