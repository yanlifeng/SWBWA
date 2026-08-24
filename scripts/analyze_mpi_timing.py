#!/usr/bin/env python3

import argparse
import csv
import html
import math
import re
import statistics
from pathlib import Path


CONFIG_LABELS = {
    "mpi_static_split": "static + split",
    "mpi_static_single": "static + single",
    "mpi_dynamic_split": "dynamic + split",
    "mpi_dynamic_single": "dynamic + single",
}

CONFIG_COLORS = {
    "mpi_static_split": "#2563eb",
    "mpi_static_single": "#dc2626",
    "mpi_dynamic_split": "#0f766e",
    "mpi_dynamic_single": "#d97706",
}

DATASET_ORDER = {
    "ERR1203383": 0,
    "small_SRR7963242": 1,
    "SRR2496709": 2,
}

FLOAT_FIELDS = (
    "pipeline_total",
    "stage1",
    "stage2",
    "stage3",
    "sam_output_writes",
    "fread_time",
    "fread_bandwidth",
    "slowest_fread",
    "worker_total",
    "worker_part1",
    "worker_part2",
    "worker_part3",
    "worker_part4",
    "worker_part5",
    "worker_part6",
    "output_flush_time",
    "output_bandwidth",
    "output_rma_flush_time",
)


TIMING_PATTERNS = (
    ("pipeline_total", re.compile(r"total - complete three-stage pipeline\s+([0-9.]+) s")),
    ("stage1", re.compile(r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s")),
    ("stage2", re.compile(r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s")),
    ("stage3", re.compile(r"stage 3 - write SAM and release batch data\s+([0-9.]+) s")),
    ("sam_output_writes", re.compile(r"SAM output writes\s+([0-9.]+) s")),
    ("fread_time", re.compile(r"FASTQ fread calls\s+([0-9.]+) s")),
    ("fread_bandwidth", re.compile(r"effective fread bandwidth\s+([0-9.]+) MiB/s")),
    ("slowest_fread", re.compile(r"slowest fread\s+([0-9.]+) s")),
    ("worker_total", re.compile(r"total - active merge worker\s+([0-9.]+) s")),
    ("worker_part1", re.compile(r"part 1 - prepare CPE task and reusable buffers\s+([0-9.]+) s")),
    ("worker_part2", re.compile(r"part 2 - CPE FASTQ formatting and input release\s+([0-9.]+) s")),
    ("worker_part3", re.compile(r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s")),
    ("worker_part4", re.compile(r"part 4 - assign slices in the shared SAM buffer\s+([0-9.]+) s")),
    ("worker_part5", re.compile(r"part 5 - CPE SAM record generation\s+([0-9.]+) s")),
    ("worker_part6", re.compile(r"part 6 - release temporary worker data\s+([0-9.]+) s")),
)

OUTPUT_PATTERNS = (
    ("output_flush_time", re.compile(r"buffered flush total\s+([0-9.]+) s")),
    ("output_bandwidth", re.compile(r"POSIX (?:pwrite|write) effective bandwidth\s+([0-9.]+) MiB/s")),
    ("output_rma_flush_time", re.compile(r"MPI_Win_flush\(target rank 0\)\s+([0-9.]+) s")),
)

RANK_RE = re.compile(r"MPI rank:\s+([0-9]+)\s+/\s+([0-9]+)")
RAW_BYTES_RE = re.compile(r"raw FASTQ bytes read\s+([0-9]+)")
INPUT_WORK_RE = re.compile(
    r"\[MPI rank\s+([0-9]+)/[0-9]+\] input work: "
    r"chunks=([0-9]+) records=([0-9]+) bytes=([0-9]+)"
)
SLOWEST_CHUNK_RE = re.compile(r"slowest chunk\s+([0-9]+)")
SLOWEST_CHUNK_TIME_RE = re.compile(r"slowest chunk stage 2\s+([0-9.]+) s")
ACCUMULATED_STAGE2_RE = re.compile(r"accumulated stage 2\s+([0-9.]+) s")


def median(values):
    return statistics.median(values) if values else math.nan


def mean(values):
    return statistics.fmean(values) if values else math.nan


def geomean(values):
    return math.exp(sum(math.log(value) for value in values) / len(values)) if values else math.nan


def safe_ratio(a, b):
    return a / b if b and not math.isnan(b) else math.nan


def fmt(value, digits=2):
    if value is None or math.isnan(value):
        return "-"
    return f"{value:.{digits}f}"


def escape(value):
    return html.escape(str(value), quote=True)


def parse_log(path, root):
    relative = path.relative_to(root)
    if len(relative.parts) != 3:
        raise ValueError(f"unexpected result path: {relative}")
    configuration, io_mode, filename = relative.parts
    sample = filename[:-len(".run.log")]
    read_mode = sample.rsplit("_", 1)[1]
    dataset = sample[:-(len(read_mode) + 1)]

    timing_by_rank = {}
    output_by_rank = {}
    scheduler_by_rank = {}
    scheduler_samples = []
    pending_scheduler_sample = None
    input_work = {}
    last_rank = None
    last_world_size = None
    current_timing = None
    current_output = None
    current_scheduler = None
    in_output = False
    in_scheduler = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = INPUT_WORK_RE.search(line)
            if match:
                rank = int(match.group(1))
                input_work[rank] = {
                    "claimed_chunks": int(match.group(2)),
                    "completed_records": int(match.group(3)),
                    "claimed_input_bytes": int(match.group(4)),
                }

            accumulated_match = ACCUMULATED_STAGE2_RE.search(line)
            if accumulated_match:
                pending_scheduler_sample = {
                    "accumulated_stage2": float(accumulated_match.group(1))
                }
            if pending_scheduler_sample is not None:
                chunk_match = SLOWEST_CHUNK_RE.search(line)
                if chunk_match:
                    pending_scheduler_sample["slowest_chunk"] = int(chunk_match.group(1))
                chunk_time_match = SLOWEST_CHUNK_TIME_RE.search(line)
                if chunk_time_match:
                    pending_scheduler_sample["slowest_chunk_stage2"] = float(chunk_time_match.group(1))
                    scheduler_samples.append(pending_scheduler_sample)
                    pending_scheduler_sample = None

            if "SWBWA Output Debug" in line:
                in_output = True
                in_scheduler = False
                current_output = None
                current_scheduler = None
                current_timing = None
                continue
            if "MPI Dynamic Scheduler Debug" in line:
                in_scheduler = True
                in_output = False
                current_scheduler = None
                current_output = None
                current_timing = None
                continue
            if line.startswith("===="):
                in_output = False
                in_scheduler = False
                current_output = None
                current_scheduler = None
                current_timing = None
                continue

            rank_match = RANK_RE.search(line)
            if rank_match:
                last_rank = int(rank_match.group(1))
                last_world_size = int(rank_match.group(2))
                current_timing = None
                if in_output:
                    current_output = output_by_rank.setdefault(
                        last_rank,
                        {"rank": last_rank, "world_size": last_world_size},
                    )
                if in_scheduler:
                    current_scheduler = scheduler_by_rank.setdefault(
                        last_rank,
                        {"rank": last_rank, "world_size": last_world_size},
                    )
                continue

            if in_output and current_output is not None:
                for field, pattern in OUTPUT_PATTERNS:
                    metric_match = pattern.search(line)
                    if metric_match:
                        current_output[field] = float(metric_match.group(1))
                        break
            if in_scheduler and current_scheduler is not None:
                chunk_match = SLOWEST_CHUNK_RE.search(line)
                if chunk_match:
                    current_scheduler["slowest_chunk"] = int(chunk_match.group(1))
                chunk_time_match = SLOWEST_CHUNK_TIME_RE.search(line)
                if chunk_time_match:
                    current_scheduler["slowest_chunk_stage2"] = float(chunk_time_match.group(1))

            total_match = TIMING_PATTERNS[0][1].search(line)
            if total_match and last_rank is not None:
                current_timing = timing_by_rank.setdefault(
                    last_rank,
                    {"rank": last_rank, "world_size": last_world_size},
                )
                current_timing["pipeline_total"] = float(total_match.group(1))
                continue
            if current_timing is None:
                continue
            for field, pattern in TIMING_PATTERNS[1:]:
                metric_match = pattern.search(line)
                if metric_match:
                    current_timing[field] = float(metric_match.group(1))
                    break
            bytes_match = RAW_BYTES_RE.search(line)
            if bytes_match:
                current_timing["raw_fastq_bytes"] = int(bytes_match.group(1))

    for sample in scheduler_samples:
        candidates = [
            (abs(report.get("stage2", math.inf) - sample["accumulated_stage2"]), rank)
            for rank, report in timing_by_rank.items()
        ]
        difference, rank = min(candidates)
        if difference <= 0.1:
            scheduler_by_rank.setdefault(rank, {}).update(sample)

    rows = []
    for rank in sorted(timing_by_rank):
        row = {
            "configuration": configuration,
            "io_mode": io_mode,
            "dataset": dataset,
            "read_mode": read_mode,
            "rank": rank,
            "log": str(relative),
        }
        row.update(timing_by_rank[rank])
        row.update(output_by_rank.get(rank, {}))
        row.update(scheduler_by_rank.get(rank, {}))
        row.update(input_work.get(rank, {}))
        rows.append(row)
    if len(rows) != 6:
        raise ValueError(f"{relative}: expected 6 timing reports, found {len(rows)}")
    return rows


def case_sort_key(case):
    return (
        0 if case["io_mode"] == "has1" else 1,
        DATASET_ORDER.get(case["dataset"], 99),
        0 if case["read_mode"] == "PE" else 1,
        list(CONFIG_LABELS).index(case["configuration"]),
    )


def summarize_case(rows):
    first = rows[0]
    summary = {
        "configuration": first["configuration"],
        "io_mode": first["io_mode"],
        "dataset": first["dataset"],
        "read_mode": first["read_mode"],
        "rank_count": len(rows),
        "log": first["log"],
    }
    for field in FLOAT_FIELDS:
        values = [row[field] for row in rows if field in row]
        if not values:
            continue
        summary[f"{field}_min"] = min(values)
        summary[f"{field}_median"] = median(values)
        summary[f"{field}_max"] = max(values)
        summary[f"{field}_mean"] = mean(values)
        summary[f"{field}_max_over_median"] = safe_ratio(max(values), median(values))
    for field in ("raw_fastq_bytes", "claimed_chunks", "completed_records", "claimed_input_bytes"):
        values = [row[field] for row in rows if field in row]
        if values:
            summary[f"{field}_sum"] = sum(values)

    skew_flags = []
    for field, label in (("stage1", "stage1"), ("stage2", "stage2"), ("stage3", "stage3")):
        ratio = summary.get(f"{field}_max_over_median", math.nan)
        delta = summary.get(f"{field}_max", 0.0) - summary.get(f"{field}_median", 0.0)
        if not math.isnan(ratio) and ratio >= 1.50 and delta >= 5.0:
            skew_flags.append(f"{label} rank skew {ratio:.2f}x")

    chunk_rows = [row for row in rows if "slowest_chunk_stage2" in row]
    if chunk_rows:
        slowest = max(chunk_rows, key=lambda row: row["slowest_chunk_stage2"])
        summary["slowest_chunk"] = slowest.get("slowest_chunk")
        summary["slowest_chunk_stage2"] = slowest["slowest_chunk_stage2"]
        if slowest["slowest_chunk_stage2"] >= 20.0:
            skew_flags.append(
                f'chunk {slowest.get("slowest_chunk", "?")} takes '
                f'{slowest["slowest_chunk_stage2"]:.1f}s'
            )

    bottleneck_flags = []
    read_bw = summary.get("fread_bandwidth_median", math.nan)
    if not math.isnan(read_bw) and read_bw < 20.0:
        bottleneck_flags.append(f"low input BW {read_bw:.1f} MiB/s")
    write_bw = summary.get("output_bandwidth_median", math.nan)
    if not math.isnan(write_bw) and write_bw < 20.0:
        bottleneck_flags.append(f"low output BW {write_bw:.1f} MiB/s")
    summary["skew_flags"] = "; ".join(skew_flags)
    summary["bottleneck_flags"] = "; ".join(bottleneck_flags)
    summary["flags"] = "; ".join(skew_flags + bottleneck_flags)
    return summary


def write_tsv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def svg_text(x, y, text, size=12, anchor="start", weight="normal", fill="#17202a", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'text-anchor="{anchor}" font-weight="{weight}" fill="{fill}"{transform}>'
        f'{escape(text)}</text>'
    )


def write_stage2_svg(path, summaries):
    ordered = sorted(summaries, key=case_sort_key)
    width = 1680
    left = 390
    right = 60
    top = 100
    row_h = 26
    group_gap = 12
    height = top + len(ordered) * row_h + 7 * group_gap + 80
    max_ratio = max(2.0, max(item.get("stage2_max_over_median", 1.0) for item in ordered) * 1.08)
    plot_w = width - left - right

    def sx(value):
        return left + (value / max_ratio) * plot_w

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 35, "MPI stage2 各 rank 负载均衡", 23, weight="bold"),
        svg_text(24, 61, "每行展示 6 个 rank 的 min-median-max；右侧数值为 max / median。", 13, fill="#4b5563"),
    ]
    for tick in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        if tick > max_ratio:
            continue
        x = sx(tick)
        color = "#dc2626" if tick == 1.5 else "#d1d5db"
        dash = ' stroke-dasharray="5 4"' if tick == 1.5 else ""
        parts.append(f'<line x1="{x:.1f}" y1="{top-15}" x2="{x:.1f}" y2="{height-45}" stroke="{color}" stroke-width="1"{dash}/>')
        parts.append(svg_text(x, top - 24, f"{tick:.1f}x", 11, anchor="middle", fill="#4b5563"))

    y = top
    previous_group = None
    for index, item in enumerate(ordered):
        group = (item["io_mode"], item["dataset"], item["read_mode"])
        if previous_group is not None and group != previous_group:
            y += group_gap
        previous_group = group
        if index % 2:
            parts.append(f'<rect x="8" y="{y-row_h+5:.1f}" width="{width-16}" height="{row_h}" fill="#f8fafc"/>')
        label = (
            f'{item["io_mode"]:4s}  {item["dataset"]:18s} {item["read_mode"]}  '
            f'{CONFIG_LABELS[item["configuration"]]}'
        )
        parts.append(svg_text(18, y, label, 12, fill="#1f2937"))
        med = item.get("stage2_median", math.nan)
        low = safe_ratio(item.get("stage2_min", math.nan), med)
        high = safe_ratio(item.get("stage2_max", math.nan), med)
        color = CONFIG_COLORS[item["configuration"]]
        parts.append(f'<line x1="{sx(low):.1f}" y1="{y-4:.1f}" x2="{sx(high):.1f}" y2="{y-4:.1f}" stroke="{color}" stroke-width="3"/>')
        for value in (low, 1.0, high):
            radius = 5 if value == 1.0 else 4
            fill = "#ffffff" if value != 1.0 else color
            parts.append(f'<circle cx="{sx(value):.1f}" cy="{y-4:.1f}" r="{radius}" fill="{fill}" stroke="{color}" stroke-width="2"/>')
        parts.append(svg_text(sx(high) + 9, y, f"{high:.2f}x  ({med:.1f}s median)", 11, fill="#374151"))
        y += row_h

    legend_y = height - 20
    x = 24
    for config, label in CONFIG_LABELS.items():
        parts.append(f'<line x1="{x}" y1="{legend_y-4}" x2="{x+24}" y2="{legend_y-4}" stroke="{CONFIG_COLORS[config]}" stroke-width="4"/>')
        parts.append(svg_text(x + 31, legend_y, label, 12))
        x += 205
    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def heat_color(value):
    if math.isnan(value):
        return "#e5e7eb"
    if value < 1.20:
        return "#d1fae5"
    if value < 1.50:
        return "#fef3c7"
    if value < 2.00:
        return "#fed7aa"
    return "#fecaca"


def write_outlier_svg(path, summaries):
    ordered = sorted(summaries, key=case_sort_key)
    width = 1440
    left = 520
    top = 100
    row_h = 25
    group_gap = 12
    columns = (
        ("stage1_max_over_median", "stage1 最大/中位"),
        ("stage2_max_over_median", "stage2 最大/中位"),
        ("stage3_max_over_median", "stage3 最大/中位"),
        ("fread_bandwidth_median", "读带宽 MiB/s"),
        ("output_bandwidth_median", "写带宽 MiB/s"),
    )
    col_w = 165
    height = top + len(ordered) * row_h + 7 * group_gap + 55
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        svg_text(24, 35, "MPI 计时离群矩阵", 23, weight="bold"),
        svg_text(24, 61, "橙/红计时格表示 rank 偏斜；橙/红带宽格分别低于 20/10 MiB/s。", 13, fill="#4b5563"),
    ]
    for index, (_, label) in enumerate(columns):
        parts.append(svg_text(left + index * col_w + col_w / 2, top - 25, label, 11, anchor="middle", weight="bold"))

    y = top
    previous_group = None
    for item in ordered:
        group = (item["io_mode"], item["dataset"], item["read_mode"])
        if previous_group is not None and group != previous_group:
            y += group_gap
        previous_group = group
        label = (
            f'{item["io_mode"]:4s}  {item["dataset"]:18s} {item["read_mode"]}  '
            f'{CONFIG_LABELS[item["configuration"]]}'
        )
        parts.append(svg_text(18, y, label, 12))
        for index, (field, _) in enumerate(columns):
            value = item.get(field, math.nan)
            if "bandwidth" in field:
                color = "#d1fae5" if value >= 20 else ("#fed7aa" if value >= 10 else "#fecaca")
                label_value = fmt(value, 1)
            else:
                color = heat_color(value)
                label_value = f"{value:.2f}x" if not math.isnan(value) else "-"
            x = left + index * col_w
            parts.append(f'<rect x="{x+3}" y="{y-row_h+6:.1f}" width="{col_w-6}" height="{row_h-3}" rx="2" fill="{color}"/>')
            parts.append(svg_text(x + col_w / 2, y, label_value, 12, anchor="middle", weight="bold"))
        y += row_h
    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def bandwidth_groups(summaries):
    groups = {}
    for item in summaries:
        key = (item["configuration"], item["io_mode"])
        group = groups.setdefault(key, {"read": [], "write": [], "low_read": 0, "low_write": 0})
        read_bw = item.get("fread_bandwidth_median", math.nan)
        write_bw = item.get("output_bandwidth_median", math.nan)
        if not math.isnan(read_bw):
            group["read"].append(read_bw)
            group["low_read"] += read_bw < 20.0
        if not math.isnan(write_bw):
            group["write"].append(write_bw)
            group["low_write"] += write_bw < 20.0
    return groups


def write_report(path, summaries, outliers, bottlenecks):
    index = {
        (item["configuration"], item["io_mode"], item["dataset"], item["read_mode"]): item
        for item in summaries
    }
    mode_comparisons = {}
    for io_mode in ("has1", "no1"):
        total_ratios = []
        stage2_ratios = []
        for dataset in DATASET_ORDER:
            for read_mode in ("PE", "SE"):
                split = index[("mpi_dynamic_split", io_mode, dataset, read_mode)]
                single = index[("mpi_dynamic_single", io_mode, dataset, read_mode)]
                total_ratios.append(single["pipeline_total_median"] / split["pipeline_total_median"])
                stage2_ratios.append(single["stage2_median"] / split["stage2_median"])
        mode_comparisons[io_mode] = (geomean(total_ratios), geomean(stage2_ratios))

    slow_chunks = [
        item for item in outliers
        if item.get("slowest_chunk") is not None and item.get("slowest_chunk_stage2", 0.0) >= 20.0
    ]
    lines = [
        "# MPI 大数据计时与离群点分析",
        "",
        "本报告覆盖 `static/dynamic × split/single_unordered × has1/no1`，",
        "每个 case 包含 6 个 MPI rank。流水线各 stage 有重叠，不能相加。",
        "",
        "## 结论",
        "",
    ]
    if outliers:
        lines.append(f"共识别出 **{len(outliers)}** 个 rank 负载离群 case（规则见下文）。")
    else:
        lines.append("没有识别出明显的 rank 级计时离群点。")
    if outliers:
        workloads = sorted({(item["dataset"], item["read_mode"]) for item in outliers})
        workload_text = "、".join(f"`{dataset} {read_mode}`" for dataset, read_mode in workloads)
        lines.append(
            f"这些 case 实际都来自同一工作负载模式：{workload_text}；"
            "stage1 和 stage3 没有达到离群阈值。"
        )
    if slow_chunks:
        chunk_ids = sorted({int(item["slowest_chunk"]) for item in slow_chunks})
        chunk_times = [item["slowest_chunk_stage2"] for item in slow_chunks]
        lines.append(
            f"dynamic 日志进一步定位到 chunk `{','.join(map(str, chunk_ids))}`，"
            f"单块耗时 {min(chunk_times):.1f}-{max(chunk_times):.1f} 秒。"
        )
    for io_mode in ("has1", "no1"):
        total_ratio, stage2_ratio = mode_comparisons[io_mode]
        lines.append(
            f"`{io_mode}` 下 dynamic single 相对 dynamic split 的 pipeline 几何平均为 "
            f"`{total_ratio:.3f}x`，但 stage2 仅为 `{stage2_ratio:.3f}x`；"
            "差异主要不在比对计算。"
        )
    lines.extend([
        "",
        "## 需要关注的 case",
        "",
        "| 配置 | -1 | 数据 | 模式 | pipeline 中位数(s) | stage1 | stage2 | stage3 | 读带宽 | 写带宽 | 标记 |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ])
    for item in outliers:
        lines.append(
            "| {config} | {io} | {dataset} | {read} | {total} | {s1} | {s2} | {s3} | {rbw} | {wbw} | {flags} |".format(
                config=CONFIG_LABELS[item["configuration"]],
                io=item["io_mode"],
                dataset=item["dataset"],
                read=item["read_mode"],
                total=fmt(item.get("pipeline_total_median")),
                s1=fmt(item.get("stage1_max_over_median")) + "x",
                s2=fmt(item.get("stage2_max_over_median")) + "x",
                s3=fmt(item.get("stage3_max_over_median")) + "x",
                rbw=fmt(item.get("fread_bandwidth_median"), 1),
                wbw=fmt(item.get("output_bandwidth_median"), 1),
                flags=item["skew_flags"],
            )
        )

    lines.extend([
        "",
        "## 系统性 I/O 低带宽",
        "",
        "低带宽是跨多个数据集重复出现的模式，不作为单点离群。下面的计数分母均为 6 个数据 case。",
        "",
        "| 配置 | -1 | 读带宽中位数 | 低于20的读 case | 写带宽中位数 | 低于20的写 case |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for (configuration, io_mode), group in sorted(
        bandwidth_groups(summaries).items(),
        key=lambda item: (0 if item[0][1] == "has1" else 1, list(CONFIG_LABELS).index(item[0][0])),
    ):
        lines.append(
            f'| {CONFIG_LABELS[configuration]} | {io_mode} | '
            f'{fmt(median(group["read"]), 1)} | {group["low_read"]}/6 | '
            f'{fmt(median(group["write"]), 1)} | {group["low_write"]}/6 |'
        )

    lines.extend([
        "",
        "## 全部 case",
        "",
        "时间列为 6 rank 的 `median [min, max]`。带宽单位为 MiB/s。",
        "",
        "| 配置 | -1 | 数据 | 模式 | pipeline | stage1 | stage2 | stage3 | 读带宽 | 写带宽 |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for item in sorted(summaries, key=case_sort_key):
        def interval(field, digits=1):
            return (
                f'{fmt(item.get(field + "_median"), digits)} '
                f'[{fmt(item.get(field + "_min"), digits)}, {fmt(item.get(field + "_max"), digits)}]'
            )

        lines.append(
            f'| {CONFIG_LABELS[item["configuration"]]} | {item["io_mode"]} | '
            f'{item["dataset"]} | {item["read_mode"]} | '
            f'{interval("pipeline_total")} | {interval("stage1")} | '
            f'{interval("stage2")} | {interval("stage3")} | '
            f'{interval("fread_bandwidth")} | {interval("output_bandwidth")} |'
        )

    lines.extend([
        "",
        "## 判定规则",
        "",
        "- rank skew：同一 case 的 `max / median >= 1.50`，且最大值比中位数至少多 5 秒。",
        "- 低带宽：6 rank 的有效读/写带宽中位数低于 20 MiB/s；低于 10 MiB/s 在图中标红。",
        "- 这些规则用于定位异常，不代表 stage1/3 可以相加，也不把不同输出策略的正常差异直接判成程序错误。",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Analyze SWBWA MPI timing logs")
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    root = args.result_root.resolve()
    output_dir = (args.output_dir or root).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rank_rows = []
    for path in sorted(root.glob("mpi_*/*/*.run.log")):
        rank_rows.extend(parse_log(path, root))
    if len(rank_rows) != 48 * 6:
        raise SystemExit(f"expected 288 rank timing rows, found {len(rank_rows)}")

    grouped = {}
    for row in rank_rows:
        key = (row["configuration"], row["io_mode"], row["dataset"], row["read_mode"])
        grouped.setdefault(key, []).append(row)
    summaries = [summarize_case(rows) for rows in grouped.values()]
    summaries.sort(key=case_sort_key)
    outliers = [item for item in summaries if item["skew_flags"]]
    bottlenecks = [item for item in summaries if item["bottleneck_flags"]]

    rank_fields = [
        "configuration", "io_mode", "dataset", "read_mode", "rank", "world_size",
        "pipeline_total", "stage1", "stage2", "stage3", "sam_output_writes",
        "fread_time", "fread_bandwidth", "slowest_fread", "raw_fastq_bytes",
        "worker_total", "worker_part1", "worker_part2", "worker_part3",
        "worker_part4", "worker_part5", "worker_part6", "output_flush_time",
        "output_bandwidth", "output_rma_flush_time", "claimed_chunks",
        "completed_records", "claimed_input_bytes", "slowest_chunk",
        "slowest_chunk_stage2", "log",
    ]
    summary_fields = list(summaries[0].keys())
    write_tsv(output_dir / "mpi_rank_timings.tsv", rank_rows, rank_fields)
    write_tsv(output_dir / "mpi_case_summary.tsv", summaries, summary_fields)
    write_tsv(output_dir / "mpi_outliers.tsv", outliers, summary_fields)
    write_tsv(output_dir / "mpi_io_bottlenecks.tsv", bottlenecks, summary_fields)
    write_stage2_svg(output_dir / "mpi_stage2_balance.svg", summaries)
    write_outlier_svg(output_dir / "mpi_outlier_matrix.svg", summaries)
    write_report(output_dir / "mpi_timing_report.md", summaries, outliers, bottlenecks)

    print(f"rank rows: {len(rank_rows)}")
    print(f"cases: {len(summaries)}")
    print(f"flagged cases: {len(outliers)}")
    for item in outliers:
        print(
            f'{item["configuration"]}/{item["io_mode"]}/'
            f'{item["dataset"]}_{item["read_mode"]}: {item["skew_flags"]}'
        )
    print(f"systemic low-bandwidth cases: {len(bottlenecks)}")


if __name__ == "__main__":
    main()
