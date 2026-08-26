#!/usr/bin/env python3
"""Summarize the CPE/LDM optimization experiment into TSV and SVG files."""

import argparse
import html
import math
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
]

PARENT = {
    "WORKER_ALIGNMENT": None,
    "MEM_CHAIN": "WORKER_ALIGNMENT",
    "MEM_CHAIN_COLLECT": "MEM_CHAIN",
    "MEM_CHAIN_BUILD": "MEM_CHAIN",
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
    "MEM_CHAIN_BUILD",
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

LADDER = [
    ("previous", "Previous optimized", "../cpe_hotspot_opt_20260825/optimized/SRR7963242_PE.log"),
    ("context", "Per-CPE context", "context_profile/SRR7963242_PE.run.log"),
    ("ksw_ldm", "KSW query profile in LDM", "ksw_ldm/SRR7963242_PE.run.log"),
    ("smem_ldm", "SMEM temporary vectors in LDM", "smem_ldm/SRR7963242_PE.run.log"),
    ("query_reuse", "Reverse-complement reuse", "ksw_query_reuse/SRR7963242_PE.run.log"),
    ("local_ksw_o3", "Local KSW DP at O3", "ksw_o3_only/SRR7963242_PE.run.log"),
    ("compact_keys", "Compact mate-dedup keys", "compact_dedup/SRR7963242_PE.run.log"),
    ("compact_keys_ldm", "Compact keys and scratch in LDM", "compact_dedup_ldm/SRR7963242_PE.run2.log"),
    ("final", "Function-scoped O3 final", "narrow_o3/SRR7963242_PE.log"),
]

BIG_CASES = [
    ("ERR1203383 PE (75 bp)", "ERR1203383_PE.run.log", "ERR1203383_PE.final.log"),
    ("small SRR7963242 PE (150 bp)", "small_SRR7963242_PE.run.log", "small_SRR7963242_PE.final.log"),
    ("small SRR7963242 SE (150 bp)", "small_SRR7963242_SE.run.log", "small_SRR7963242_SE.log"),
    ("SRR2496709 PE", "SRR2496709_PE.run.log", "SRR2496709_PE.final.log"),
]

SCALED_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]*)?)([KMGT]?)\s*$")


def scaled_number(value):
    match = SCALED_RE.match(value)
    if not match:
        raise ValueError("invalid LWPF number: {!r}".format(value))
    scale = {"": 1.0, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    return float(match.group(1)) * scale[match.group(2)]


def last_float(pattern, text):
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if not matches:
        raise ValueError("missing timing pattern: {}".format(pattern))
    return float(matches[-1])


def parse_timing(path):
    text = path.read_text(errors="replace")
    return {
        "stage1": last_float(
            r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s", text),
        "stage2": last_float(
            r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s", text),
        "stage3": last_float(
            r"stage 3 - write SAM and release batch data\s+([0-9.]+) s", text),
        "part3": last_float(
            r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s", text),
        "read_bw": last_float(
            r"effective fread bandwidth\s+([0-9.]+) MiB/s", text),
    }


def parse_lwpf(path):
    text = path.read_text(errors="replace")
    rows = []
    in_table = False
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\r")
        if line.startswith("LWPF kernel summary"):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("===="):
            break
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.split("|")[1:-1]]
        if len(cells) != 25:
            continue
        values = [scaled_number(cell) for cell in cells[1:]]
        rows.append(values)
    if len(rows) != len(PROFILE_REGIONS):
        raise ValueError("{}: expected {} LWPF rows, found {}".format(
            path, len(PROFILE_REGIONS), len(rows)))
    result = {}
    for name, values in zip(PROFILE_REGIONS, rows):
        result[name] = {
            "cycles_avg": values[0],
            "cycles_min": values[1],
            "cycles_max": values[2],
            "instructions_avg": values[3],
            "dcache_accesses_avg": values[12],
            "dcache_misses_avg": values[15],
            "instruction_stall_cycles_avg": values[21],
        }
    return result


def tree_depth(region):
    value = 0
    parent = PARENT[region]
    while parent is not None:
        value += 1
        parent = PARENT[parent]
    return value


def svg_start(width, height, title, subtitle):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
        'viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text { font-family: -apple-system, BlinkMacSystemFont, '
        '"Segoe UI", Arial, sans-serif; fill: #17242c; letter-spacing: 0; } '
        '.small { font-size: 12px; fill: #59666e; } '
        '.mono { font-family: Menlo, Consolas, monospace; }</style>',
        '<text x="36" y="39" font-size="22" font-weight="700">{}</text>'.format(
            html.escape(title)),
        '<text x="36" y="62" class="small">{}</text>'.format(
            html.escape(subtitle)),
    ]


def write_ladder_tsv(rows, path):
    lines = ["variant\tlabel\tstage2_s\tpart3_s\tstage2_speedup_vs_previous\tpart3_speedup_vs_previous"]
    first = rows[0][2]
    for variant, label, timing in rows:
        lines.append("\t".join([
            variant, label, "{:.3f}".format(timing["stage2"]),
            "{:.3f}".format(timing["part3"]),
            "{:.4f}".format(first["stage2"] / timing["stage2"]),
            "{:.4f}".format(first["part3"] / timing["part3"]),
        ]))
    path.write_text("\n".join(lines) + "\n")


def write_ladder_svg(rows, path):
    width, row_height = 1220, 45
    height = 108 + len(rows) * row_height
    left, bar_width = 320, 720
    maximum = max(item[2]["stage2"] for item in rows)
    out = svg_start(
        width, height, "CPE optimization ladder: SRR7963242 PE subset",
        "One CGS-cross process; lower is better; Stage 2 includes part 3",
    )
    out.extend([
        '<rect x="865" y="28" width="13" height="9" fill="#697982"/>',
        '<text x="884" y="37" class="small">Stage 2</text>',
        '<rect x="950" y="28" width="13" height="9" fill="#d55e00"/>',
        '<text x="969" y="37" class="small">CPE part 3</text>',
    ])
    baseline = rows[0][2]["stage2"]
    for index, (_, label, timing) in enumerate(rows):
        y = 82 + index * row_height
        stage2_width = bar_width * timing["stage2"] / maximum
        part3_width = bar_width * timing["part3"] / maximum
        gain = (1.0 - timing["stage2"] / baseline) * 100.0
        out.append('<text x="36" y="{}" font-size="12">{}</text>'.format(
            y + 16, html.escape(label)))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="9" fill="#697982"/>'.format(
            left, y, stage2_width))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="9" fill="#d55e00"/>'.format(
            left, y + 12, part3_width))
        out.append('<text x="{}" y="{}" class="small">{:.3f} s ({:+.1f}%)</text>'.format(
            left + stage2_width + 8, y + 9, timing["stage2"], -gain))
    out.append("</svg>")
    path.write_text("\n".join(out) + "\n")


def write_bigdata_tsv(rows, path):
    lines = [
        "case\tbaseline_stage2_s\tfinal_stage2_s\tstage2_speedup\t"
        "baseline_part3_s\tfinal_part3_s\tpart3_speedup\tfinal_read_bw_mib_s"
    ]
    for label, baseline, final in rows:
        lines.append("\t".join([
            label, "{:.3f}".format(baseline["stage2"]),
            "{:.3f}".format(final["stage2"]),
            "{:.4f}".format(baseline["stage2"] / final["stage2"]),
            "{:.3f}".format(baseline["part3"]),
            "{:.3f}".format(final["part3"]),
            "{:.4f}".format(baseline["part3"] / final["part3"]),
            "{:.3f}".format(final["read_bw"]),
        ]))
    path.write_text("\n".join(lines) + "\n")


def write_bigdata_svg(rows, path):
    width, row_height = 1260, 58
    height = 112 + len(rows) * row_height
    left, bar_width = 340, 720
    maximum = max(item[1]["stage2"] for item in rows)
    out = svg_start(
        width, height, "Big-data CPE alignment: previous versus final",
        "Stage 2 and part 3; final compute-only runs write to /dev/null",
    )
    out.extend([
        '<rect x="860" y="28" width="13" height="9" fill="#697982"/>',
        '<text x="879" y="37" class="small">previous Stage 2</text>',
        '<rect x="1000" y="28" width="13" height="9" fill="#d55e00"/>',
        '<text x="1019" y="37" class="small">final Stage 2</text>',
    ])
    for index, (label, baseline, final) in enumerate(rows):
        y = 83 + index * row_height
        old_width = bar_width * baseline["stage2"] / maximum
        new_width = bar_width * final["stage2"] / maximum
        speedup = baseline["stage2"] / final["stage2"]
        out.append('<text x="36" y="{}" font-size="12">{}</text>'.format(
            y + 16, html.escape(label)))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="10" fill="#697982"/>'.format(
            left, y, old_width))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="10" fill="#d55e00"/>'.format(
            left, y + 14, new_width))
        out.append('<text x="{}" y="{}" class="small">{:.3f} s, {:.3f}x; part3 {:.3f} s</text>'.format(
            left + new_width + 8, y + 23, final["stage2"], speedup,
            final["part3"]))
    out.append("</svg>")
    path.write_text("\n".join(out) + "\n")


def write_hotspot_tsv(profiles, path):
    lines = [
        "profile\tregion\tparent\tdepth\tcycles_avg\tcycles_min\tcycles_max\t"
        "worker_cycle_percent\tinstructions_avg\tdcache_accesses_avg\t"
        "dcache_misses_avg\tdcache_miss_percent\tinstruction_stall_cycles_avg"
    ]
    for profile_name, regions in profiles.items():
        worker_cycles = regions["WORKER_ALIGNMENT"]["cycles_avg"]
        for region in TREE_ORDER:
            values = regions[region]
            accesses = values["dcache_accesses_avg"]
            miss_rate = (100.0 * values["dcache_misses_avg"] / accesses
                         if accesses else math.nan)
            lines.append("\t".join([
                profile_name, region, PARENT[region] or "", str(tree_depth(region)),
                "{:.0f}".format(values["cycles_avg"]),
                "{:.0f}".format(values["cycles_min"]),
                "{:.0f}".format(values["cycles_max"]),
                "{:.4f}".format(100.0 * values["cycles_avg"] / worker_cycles),
                "{:.0f}".format(values["instructions_avg"]),
                "{:.0f}".format(accesses),
                "{:.0f}".format(values["dcache_misses_avg"]),
                "{:.6f}".format(miss_rate),
                "{:.0f}".format(values["instruction_stall_cycles_avg"]),
            ]))
    path.write_text("\n".join(lines) + "\n")


def write_hotspot_svg(profiles, path):
    width, row_height = 1510, 29
    height = 105 + len(TREE_ORDER) * row_height
    label_x, bar_x, bar_width = 36, 390, 880
    colors = {"ERR 75bp PE": "#0072b2", "SRR 150bp PE": "#d55e00"}
    out = svg_start(
        width, height, "Final CPE kernel hierarchy by read length",
        "Bars are percentages of WORKER_ALIGNMENT; indented child regions are contained in parents and must not be added",
    )
    legend_x = 1030
    for index, name in enumerate(profiles):
        x = legend_x + index * 150
        out.append('<rect x="{}" y="28" width="13" height="9" fill="{}"/>'.format(
            x, colors[name]))
        out.append('<text x="{}" y="37" class="small">{}</text>'.format(
            x + 19, html.escape(name)))
    for index, region in enumerate(TREE_ORDER):
        y = 79 + index * row_height
        level = tree_depth(region)
        label = region.replace("_", " ").title()
        if level:
            out.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#c8ced3"/>'.format(
                label_x + level * 14 - 7, y - 9, label_x + level * 14 - 7, y + 8))
        out.append('<text x="{}" y="{}" font-size="12">{}</text>'.format(
            label_x + level * 14, y + 8, html.escape(label)))
        for profile_index, (name, regions) in enumerate(profiles.items()):
            worker = regions["WORKER_ALIGNMENT"]["cycles_avg"]
            percent = 100.0 * regions[region]["cycles_avg"] / worker
            length = bar_width * percent / 100.0
            bar_y = y - 4 + profile_index * 9
            out.append('<rect x="{}" y="{}" width="{:.2f}" height="7" fill="{}"/>'.format(
                bar_x, bar_y, length, colors[name]))
            out.append('<text x="{}" y="{}" class="small">{:.1f}%</text>'.format(
                bar_x + length + 6, bar_y + 7, percent))
    out.append("</svg>")
    path.write_text("\n".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path,
        default=Path("correctness_results/cpe_ldm_opt_20260826"))
    parser.add_argument(
        "--baseline-root", type=Path,
        default=Path("correctness_results/cpe_bigdata_opt_20260825/has1"))
    args = parser.parse_args()

    root = args.root
    output = root / "analysis"
    output.mkdir(parents=True, exist_ok=True)

    ladder = []
    for variant, label, relative in LADDER:
        path = root / relative
        ladder.append((variant, label, parse_timing(path)))
    write_ladder_tsv(ladder, output / "optimization_ladder.tsv")
    write_ladder_svg(ladder, output / "optimization_ladder.svg")

    bigdata = []
    for label, baseline_name, final_name in BIG_CASES:
        baseline = parse_timing(args.baseline_root / baseline_name)
        final = parse_timing(root / "narrow_o3" / final_name)
        bigdata.append((label, baseline, final))
    write_bigdata_tsv(bigdata, output / "bigdata_stage2.tsv")
    write_bigdata_svg(bigdata, output / "bigdata_stage2.svg")

    profiles = {
        "ERR 75bp PE": parse_lwpf(root / "narrow_o3" / "ERR1203383_PE.profile.log"),
        "SRR 150bp PE": parse_lwpf(root / "narrow_o3" / "SRR7963242_PE.profile.log"),
    }
    write_hotspot_tsv(profiles, output / "final_hotspots.tsv")
    write_hotspot_svg(profiles, output / "final_hotspots.svg")


if __name__ == "__main__":
    main()
