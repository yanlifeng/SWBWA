#!/usr/bin/env python3

import csv
import math
import os
import re
import statistics


ROOT = os.path.dirname(os.path.abspath(__file__))
LINE_RE = re.compile(
    r"Processed (\d+) reads in ([0-9.]+) CPU sec, ([0-9.]+) real sec"
)

DATASETS = [
    {
        "name": "ERR1203383",
        "color": "#c43c39",
        "k": 716400,
        "run1": "swbwa_profile_ERR1203383_1MiB.log",
        "run2": "swbwa_profile_ERR1203383_1MiB_repeat.log",
    },
    {
        "name": "SRR2496709",
        "color": "#2563a6",
        "k": 772300,
        "run1": "swbwa_profile_SRR2496709_1MiB.log",
        "run2": "swbwa_profile_SRR2496709_1MiB_repeat.log",
    },
    {
        "name": "small_SRR7963242",
        "color": "#23845f",
        "k": 846600,
        "run1": "swbwa_profile_small_SRR7963242_1MiB.log",
        "run2": None,
    },
]


def parse_log(filename):
    values = []
    with open(os.path.join(ROOT, filename), encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = LINE_RE.search(line)
            if match:
                values.append(
                    {
                        "reads": int(match.group(1)),
                        "cpu": float(match.group(2)),
                        "real": float(match.group(3)),
                    }
                )

    full_reads = max(item["reads"] for item in values)
    for item in values:
        scale = full_reads / item["reads"]
        item["cpu_normalized"] = item["cpu"] * scale
        item["real_normalized"] = item["real"] * scale
    return values, full_reads


def quantile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    fraction = position - lower
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def rolling_median(values, width):
    radius = width // 2
    return [
        statistics.median(values[max(0, i - radius):min(len(values), i + radius + 1)])
        for i in range(len(values))
    ]


def xml_escape(value):
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def svg_text(parts, x, y, value, size=16, fill="#1f2933", anchor="start", weight="400"):
    parts.append(
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
        f'font-weight="{weight}">{xml_escape(value)}</text>'
    )


def make_polyline(values, x0, y0, width, height, y_min, y_max):
    log_min = math.log10(y_min)
    log_span = math.log10(y_max) - log_min
    count = max(1, len(values) - 1)
    points = []
    for i, value in enumerate(values):
        x = x0 + width * i / count
        clipped = min(y_max, max(y_min, value))
        y = y0 + height * (1.0 - (math.log10(clipped) - log_min) / log_span)
        points.append(f"{x:.2f},{y:.2f}")
    return " ".join(points)


def render_svg(results):
    canvas_width = 1600
    canvas_height = 1120
    plot_x = 125
    plot_width = 1370
    panel_height = 250
    panel_gap = 70
    first_y = 175
    y_min = 0.25
    y_max = 40.0
    y_ticks = [0.5, 1, 2, 5, 10, 20, 40]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas_width}" '
        f'height="{canvas_height}" viewBox="0 0 {canvas_width} {canvas_height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]

    svg_text(parts, 70, 55, "BWA-MEM compute cost by approximately 1 MiB FASTQ chunk", 30, weight="700")
    svg_text(
        parts,
        70,
        88,
        "PE input; 64 threads; CPU seconds normalized to a full chunk; logarithmic y-axis",
        17,
        fill="#52606d",
    )
    parts.append('<line x1="75" y1="118" x2="125" y2="118" stroke="#64748b" stroke-width="1.2" opacity="0.35"/>')
    svg_text(parts, 135, 124, "individual run", 15, fill="#52606d")
    parts.append('<line x1="290" y1="118" x2="350" y2="118" stroke="#111827" stroke-width="3"/>')
    svg_text(parts, 360, 124, "21-chunk rolling median of reproducible mean", 15, fill="#52606d")

    for row, result in enumerate(results):
        y0 = first_y + row * (panel_height + panel_gap)
        values = result["mean_cpu"]
        smooth = rolling_median(values, 21)
        color = result["color"]

        parts.append(
            f'<rect x="{plot_x}" y="{y0}" width="{plot_width}" height="{panel_height}" '
            'fill="#fbfcfd" stroke="#cbd5e1" stroke-width="1"/>'
        )
        if result["name"] == "ERR1203383":
            tail_x = plot_x + plot_width * 0.99
            parts.append(
                f'<rect x="{tail_x:.2f}" y="{y0}" width="{plot_x + plot_width - tail_x:.2f}" '
                f'height="{panel_height}" fill="#fee2e2" opacity="0.75"/>'
            )

        for tick in y_ticks:
            log_min = math.log10(y_min)
            log_span = math.log10(y_max) - log_min
            tick_y = y0 + panel_height * (
                1.0 - (math.log10(tick) - log_min) / log_span
            )
            parts.append(
                f'<line x1="{plot_x}" y1="{tick_y:.2f}" x2="{plot_x + plot_width}" '
                f'y2="{tick_y:.2f}" stroke="#dbe3ea" stroke-width="1"/>'
            )
            svg_text(parts, plot_x - 14, tick_y + 5, f"{tick:g}", 14, fill="#52606d", anchor="end")

        for tick in range(0, 101, 20):
            tick_x = plot_x + plot_width * tick / 100
            parts.append(
                f'<line x1="{tick_x:.2f}" y1="{y0}" x2="{tick_x:.2f}" '
                f'y2="{y0 + panel_height}" stroke="#e5eaf0" stroke-width="1"/>'
            )
            svg_text(parts, tick_x, y0 + panel_height + 25, f"{tick}%", 14, fill="#52606d", anchor="middle")

        for run in result["cpu_runs"]:
            points = make_polyline(run, plot_x, y0, plot_width, panel_height, y_min, y_max)
            parts.append(
                f'<polyline points="{points}" fill="none" stroke="{color}" '
                'stroke-width="0.8" opacity="0.28" stroke-linejoin="round"/>'
            )

        smooth_points = make_polyline(smooth, plot_x, y0, plot_width, panel_height, y_min, y_max)
        parts.append(
            f'<polyline points="{smooth_points}" fill="none" stroke="{color}" '
            'stroke-width="3" stroke-linejoin="round"/>'
        )

        stats = result["cpu_stats"]
        title = (
            f'{result["name"]}  |  {result["chunks"]} chunks  |  '
            f'{result["full_reads"] // 2} pairs/chunk  |  -K {result["k"]}'
        )
        svg_text(parts, plot_x + 12, y0 + 24, title, 17, color, weight="700")
        summary = (
            f'median {stats["median"]:.3f}s   p95 {stats["p95"]:.3f}s   '
            f'p99 {stats["p99"]:.3f}s   max {stats["max"]:.3f}s'
        )
        svg_text(parts, plot_x + plot_width - 12, y0 + 24, summary, 15, fill="#334e68", anchor="end")

    svg_text(parts, canvas_width / 2, canvas_height - 30, "Position in each FASTQ file", 17, fill="#334e68", anchor="middle", weight="700")
    parts.append('</svg>')
    with open(os.path.join(ROOT, "bwamem_1mib_cpu_profile.svg"), "w", encoding="utf-8") as stream:
        stream.write("\n".join(parts))


def main():
    results = []
    csv_rows = []

    for config in DATASETS:
        run1, full_reads = parse_log(config["run1"])
        runs = [run1]
        if config["run2"]:
            run2, repeat_full_reads = parse_log(config["run2"])
            if len(run2) != len(run1) or repeat_full_reads != full_reads:
                raise RuntimeError(f'incompatible repeat for {config["name"]}')
            runs.append(run2)

        cpu_runs = [[item["cpu_normalized"] for item in run] for run in runs]
        real_runs = [[item["real_normalized"] for item in run] for run in runs]
        mean_cpu = [statistics.mean(values) for values in zip(*cpu_runs)]
        mean_real = [statistics.mean(values) for values in zip(*real_runs)]
        cpu_stats = {
            "mean": statistics.mean(mean_cpu),
            "median": quantile(mean_cpu, 0.50),
            "p95": quantile(mean_cpu, 0.95),
            "p99": quantile(mean_cpu, 0.99),
            "max": max(mean_cpu),
        }
        real_stats = {
            "mean": statistics.mean(mean_real),
            "median": quantile(mean_real, 0.50),
            "p95": quantile(mean_real, 0.95),
            "p99": quantile(mean_real, 0.99),
            "max": max(mean_real),
        }
        result = dict(config)
        result.update(
            {
                "chunks": len(run1),
                "full_reads": full_reads,
                "cpu_runs": cpu_runs,
                "real_runs": real_runs,
                "mean_cpu": mean_cpu,
                "mean_real": mean_real,
                "cpu_stats": cpu_stats,
                "real_stats": real_stats,
            }
        )
        results.append(result)

        for index in range(len(run1)):
            row = {
                "dataset": config["name"],
                "chunk_index": index,
                "file_position_percent": 100.0 * index / max(1, len(run1) - 1),
                "read_ends": run1[index]["reads"],
                "cpu_seconds_run1": cpu_runs[0][index],
                "real_seconds_run1": real_runs[0][index],
                "cpu_seconds_run2": cpu_runs[1][index] if len(runs) > 1 else "",
                "real_seconds_run2": real_runs[1][index] if len(runs) > 1 else "",
                "cpu_seconds_mean": mean_cpu[index],
                "real_seconds_mean": mean_real[index],
            }
            csv_rows.append(row)

    csv_path = os.path.join(ROOT, "bwamem_1mib_chunk_profile.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)

    render_svg(results)

    print("dataset\tchunks\tcpu_median\tcpu_p95\tcpu_p99\tcpu_max\treal_median\treal_p95\treal_p99\treal_max")
    for result in results:
        cpu = result["cpu_stats"]
        real = result["real_stats"]
        print(
            f'{result["name"]}\t{result["chunks"]}\t{cpu["median"]:.4f}\t'
            f'{cpu["p95"]:.4f}\t{cpu["p99"]:.4f}\t{cpu["max"]:.4f}\t'
            f'{real["median"]:.4f}\t{real["p95"]:.4f}\t{real["p99"]:.4f}\t{real["max"]:.4f}'
        )


if __name__ == "__main__":
    main()
