#!/usr/bin/env python3
"""Parse SWBWA LWPF logs and generate reproducible TSV/SVG reports."""

import argparse
import html
import math
import re
from pathlib import Path


REGIONS = [
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
]

METRICS = [
    "cycles",
    "instructions",
    "global_loads",
    "global_stores",
    "dcache_accesses",
    "dcache_misses",
    "memory_barrier_wait_cycles",
    "instruction_buffer_empty_cycles",
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
    "PAIRING",
    "SAM_COPY",
]

STAGE2_PARTS = [
    "prepare CPE task and reusable buffers",
    "CPE FASTQ formatting and input release",
    "CPE alignment and SAM length pass",
    "assign slices in the shared SAM buffer",
    "CPE SAM record generation",
    "release temporary worker data",
]

EXPECTED_MD5 = {
    ("ERR1203383", "PE"): "c2af4bf0b057d5125ce2a9d770f13741",
    ("ERR1203383", "SE"): "473eec3972fbd651d8b911b9fb5c6e25",
    ("SRR2496709", "PE"): "20f980ce3955d09e7c133ba26ddf2b77",
    ("SRR7963242", "PE"): "a799dc7268f389120ec92820e04b5118",
    ("SRR7963242", "SE"): "0acfb1f46abd9fed3862c28a33e26da4",
}

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
NUMBER_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]*)?)([KMGT]?)\s*$")


def parse_scaled_number(text):
    match = NUMBER_RE.match(ANSI_RE.sub("", text))
    if not match:
        raise ValueError("invalid LWPF value: {!r}".format(text))
    scale = {"": 1.0, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
    return float(match.group(1)) * scale[match.group(2)]


def parse_case_name(path):
    match = re.match(r"(.+)_(PE|SE)\.log$", path.name)
    if not match:
        return None
    return match.group(1), match.group(2)


def parse_lwpf(text):
    rows = []
    in_table = False
    for raw_line in text.splitlines():
        line = ANSI_RE.sub("", raw_line.rstrip("\r"))
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
        values = [parse_scaled_number(cell) for cell in cells[1:]]
        rows.append(values)

    if not rows:
        raise ValueError("LWPF table not found")
    if len(rows) not in (len(REGIONS) - 1, len(REGIONS)):
        raise ValueError("unexpected LWPF region count: {}".format(len(rows)))

    names = REGIONS[:len(rows)]
    result = {}
    for region, values in zip(names, rows):
        metrics = {}
        for index, metric in enumerate(METRICS):
            offset = index * 3
            metrics[metric] = {
                "avg": values[offset],
                "min": values[offset + 1],
                "max": values[offset + 2],
            }
        result[region] = metrics
    return result


def parse_timing(text):
    result = {"pipeline": {}, "stage2_parts": {}}
    for line in text.splitlines():
        match = re.match(r"\s*stage ([123]) - .*?\s+([0-9.]+) s\s*$", line)
        if match:
            result["pipeline"]["stage{}".format(match.group(1))] = float(match.group(2))

    marker = text.find("Stage 2 worker details")
    if marker < 0:
        return result
    section = text[marker:]
    match = re.search(r"total - active merge worker\s+([0-9.]+) s", section)
    if match:
        result["stage2_total"] = float(match.group(1))
    for index, label in enumerate(STAGE2_PARTS, 1):
        pattern = r"part {} - {}\s+([0-9.]+) s".format(index, re.escape(label))
        match = re.search(pattern, section)
        if match:
            result["stage2_parts"]["part{}".format(index)] = float(match.group(1))
    return result


def depth(region):
    value = 0
    parent = PARENT[region]
    while parent is not None:
        value += 1
        parent = PARENT[parent]
    return value


def discover_logs(root):
    records = []
    for path in sorted(root.glob("*/*.log")):
        case = parse_case_name(path)
        if case is None:
            continue
        text = path.read_text(errors="replace")
        dataset, read_mode = case
        records.append({
            "variant": path.parent.name,
            "dataset": dataset,
            "read_mode": read_mode,
            "path": path,
            "lwpf": parse_lwpf(text),
            "timing": parse_timing(text),
        })
    return records


def write_kernel_tsv(records, path):
    columns = ["variant", "dataset", "read_mode", "region", "parent", "depth"]
    for metric in METRICS:
        columns.extend([metric + "_avg", metric + "_min", metric + "_max"])
    columns.append("cpi_avg")
    lines = ["\t".join(columns)]
    for record in records:
        for region in REGIONS:
            if region not in record["lwpf"]:
                continue
            row = [
                record["variant"], record["dataset"], record["read_mode"],
                region, PARENT[region] or "", str(depth(region)),
            ]
            for metric in METRICS:
                values = record["lwpf"][region][metric]
                row.extend("{:.0f}".format(values[key]) for key in ("avg", "min", "max"))
            cycles = record["lwpf"][region]["cycles"]["avg"]
            instructions = record["lwpf"][region]["instructions"]["avg"]
            row.append("{:.6f}".format(cycles / instructions)
                       if instructions else "nan")
            lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")


def write_stage2_tsv(records, path):
    columns = [
        "variant", "dataset", "read_mode", "stage1_seconds",
        "stage2_seconds", "stage3_seconds", "worker_total_seconds",
    ] + ["part{}_seconds".format(index) for index in range(1, 7)]
    lines = ["\t".join(columns)]
    for record in records:
        timing = record["timing"]
        pipeline = timing["pipeline"]
        parts = timing["stage2_parts"]
        values = [
            record["variant"], record["dataset"], record["read_mode"],
            "{:.6f}".format(pipeline.get("stage1", math.nan)),
            "{:.6f}".format(pipeline.get("stage2", math.nan)),
            "{:.6f}".format(pipeline.get("stage3", math.nan)),
            "{:.6f}".format(timing.get("stage2_total", math.nan)),
        ]
        values.extend("{:.6f}".format(parts.get("part{}".format(index), math.nan))
                      for index in range(1, 7))
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n")


def write_correctness_tsv(root, path):
    lines = ["variant\tdataset\tread_mode\tobserved_md5\texpected_md5\tstatus"]
    for md5_path in sorted(root.glob("*/*.md5")):
        match = re.match(r"(.+)_(PE|SE)\.md5$", md5_path.name)
        if not match:
            continue
        dataset, read_mode = match.groups()
        observed = md5_path.read_text().split()[0]
        expected = EXPECTED_MD5.get((dataset, read_mode), "")
        status = "PASS" if expected and observed == expected else "UNKNOWN"
        lines.append("\t".join([
            md5_path.parent.name, dataset, read_mode, observed, expected, status,
        ]))
    path.write_text("\n".join(lines) + "\n")


def selected_record(records, dataset, variant):
    candidates = [record for record in records
                  if record["dataset"] == dataset and record["read_mode"] == "PE"
                  and record["variant"] == variant]
    return candidates[0] if candidates else None


def baseline_record(records, dataset):
    return (selected_record(records, dataset, "baseline_dedup") or
            selected_record(records, dataset, "baseline"))


def svg_header(width, height, title):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
        'viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", '
        'Arial, sans-serif; fill: #1d252c; } .mono { font-family: Menlo, Consolas, '
        'monospace; } .small { font-size: 12px; fill: #53606a; }</style>',
        '<text x="40" y="42" font-size="24" font-weight="700">{}</text>'.format(
            html.escape(title)),
    ]


def write_hierarchy_svg(records, path):
    datasets = ["ERR1203383", "SRR7963242"]
    width = 1500
    left = 430
    bar_width = 950
    row_height = 27
    panel_height = 70 + len(TREE_ORDER) * row_height
    height = 90 + panel_height * len(datasets) + 60
    out = svg_header(width, height, "CPE kernel hierarchy: average cycles per sampled CPE")
    out.append('<rect x="1010" y="21" width="14" height="10" fill="#657786"/>')
    out.append('<text x="1032" y="31" class="small">baseline</text>')
    out.append('<rect x="1130" y="21" width="14" height="10" fill="#d55e00"/>')
    out.append('<text x="1152" y="31" class="small">optimized</text>')

    for panel_index, dataset in enumerate(datasets):
        baseline = baseline_record(records, dataset)
        optimized = selected_record(records, dataset, "optimized")
        if baseline is None or optimized is None:
            continue
        top = 70 + panel_index * panel_height
        max_cycles = max(
            baseline["lwpf"]["WORKER_ALIGNMENT"]["cycles"]["avg"],
            optimized["lwpf"]["WORKER_ALIGNMENT"]["cycles"]["avg"],
        )
        out.append('<text x="40" y="{}" font-size="18" font-weight="700">{} PE</text>'.format(
            top + 20, dataset))
        out.append('<text x="{}" y="{}" class="small">0</text>'.format(left, top + 20))
        out.append('<text x="{}" y="{}" class="small" text-anchor="end">{:.1f}G cycles</text>'.format(
            left + bar_width, top + 20, max_cycles / 1e9))
        out.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#c8ced3"/>'.format(
            left, top + 28, left + bar_width, top + 28))

        row_y = {}
        for row_index, region in enumerate(TREE_ORDER):
            y = top + 55 + row_index * row_height
            row_y[region] = y
            level = depth(region)
            label_x = 50 + level * 23
            has_children = any(parent == region for parent in PARENT.values())
            if PARENT[region] is not None:
                parent_y = row_y[PARENT[region]]
                connector_x = label_x - 13
                out.append('<path d="M {} {} V {} H {}" fill="none" stroke="#aab3ba" '
                           'stroke-width="1"/>'.format(connector_x, parent_y + 7,
                                                       y - 3, label_x - 4))
            weight = "700" if has_children or PARENT[region] is None else "400"
            out.append('<text x="{}" y="{}" class="mono" font-size="12" '
                       'font-weight="{}">{}</text>'.format(
                           label_x, y + 4, weight, html.escape(region)))

            baseline_has_region = region in baseline["lwpf"]
            baseline_cycles = baseline["lwpf"].get(region, {}).get("cycles", {}).get("avg", 0.0)
            optimized_cycles = optimized["lwpf"].get(region, {}).get("cycles", {}).get("avg", 0.0)
            baseline_length = bar_width * baseline_cycles / max_cycles
            optimized_length = bar_width * optimized_cycles / max_cycles
            if baseline_has_region:
                out.append('<rect x="{}" y="{}" width="{:.2f}" height="8" fill="#657786"/>'.format(
                    left, y - 8, baseline_length))
            else:
                out.append('<text x="{}" y="{}" class="small">baseline n/a</text>'.format(
                    left + 5, y - 1))
            out.append('<rect x="{}" y="{}" width="{:.2f}" height="8" fill="#d55e00"/>'.format(
                left, y + 2, optimized_length))
            value = optimized_cycles / 1e9
            out.append('<text x="{}" y="{}" class="small">{:.2f}G</text>'.format(
                min(left + optimized_length + 6, width - 75), y + 10, value))

    out.append('<text x="40" y="{}" class="small">Nested regions overlap with their parents; '
               'do not add parent and child cycle counts.</text>'.format(height - 25))
    out.append('</svg>')
    path.write_text("\n".join(out) + "\n")


def write_stage2_svg(records, path):
    datasets = ["ERR1203383", "SRR7963242"]
    width, height = 1450, 500
    left, chart_width = 260, 1080
    colors = ["#e69f00", "#009e73", "#0072b2", "#d55e00", "#8b5fbf", "#7a7a7a"]
    out = svg_header(width, height, "Stage2 breakdown: baseline versus optimized")
    legend_x = 280
    for index, (label, color) in enumerate(zip(STAGE2_PARTS, colors)):
        x = legend_x + (index % 3) * 360
        y = 70 + (index // 3) * 24
        out.append('<rect x="{}" y="{}" width="12" height="12" fill="{}"/>'.format(x, y - 10, color))
        out.append('<text x="{}" y="{}" class="small">part{} {}</text>'.format(
            x + 18, y, index + 1, html.escape(label)))

    comparisons = []
    for dataset in datasets:
        baseline = baseline_record(records, dataset)
        optimized = selected_record(records, dataset, "optimized")
        if baseline and optimized:
            comparisons.extend([(dataset + " baseline", baseline),
                                (dataset + " optimized", optimized)])
    max_total = max(record["timing"].get("stage2_total", 0.0)
                    for _, record in comparisons)
    for row_index, (label, record) in enumerate(comparisons):
        y = 150 + row_index * 75
        out.append('<text x="40" y="{}" font-size="14">{}</text>'.format(
            y + 16, html.escape(label)))
        x = left
        for index, color in enumerate(colors, 1):
            seconds = record["timing"]["stage2_parts"].get("part{}".format(index), 0.0)
            segment = chart_width * seconds / max_total
            out.append('<rect x="{:.2f}" y="{}" width="{:.2f}" height="28" fill="{}"/>'.format(
                x, y, segment, color))
            x += segment
        total = record["timing"].get("stage2_total", 0.0)
        out.append('<text x="{}" y="{}" font-size="13" font-weight="700">{:.3f} s</text>'.format(
            min(x + 8, width - 80), y + 19, total))
    out.append('<line x1="{}" y1="130" x2="{}" y2="130" stroke="#c8ced3"/>'.format(
        left, left + chart_width))
    out.append('<text x="{}" y="125" class="small">0</text>'.format(left))
    out.append('<text x="{}" y="125" class="small" text-anchor="end">{:.1f} s</text>'.format(
        left + chart_width, max_total))
    out.append('</svg>')
    path.write_text("\n".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("correctness_results/cpe_hotspot_opt_20260825"))
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or args.root / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    records = discover_logs(args.root)
    if not records:
        raise SystemExit("no profiling logs found under {}".format(args.root))
    write_kernel_tsv(records, output_dir / "kernel_profile.tsv")
    write_stage2_tsv(records, output_dir / "stage2_profile.tsv")
    write_correctness_tsv(args.root, output_dir / "correctness.tsv")
    write_hierarchy_svg(records, output_dir / "cpe_kernel_hierarchy.svg")
    write_stage2_svg(records, output_dir / "stage2_breakdown.svg")
    print("parsed {} logs; results written to {}".format(len(records), output_dir))


if __name__ == "__main__":
    main()
