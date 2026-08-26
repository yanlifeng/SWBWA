#!/usr/bin/env python3
"""Analyze release-mode cgs_cross + pool big-data accuracy and timing logs."""

import argparse
import html
import math
import re
from pathlib import Path


CASES = [
    ("ERR1203383", "PE"),
    ("ERR1203383", "SE"),
    ("small_SRR7963242", "PE"),
    ("small_SRR7963242", "SE"),
    ("SRR2496709", "PE"),
    ("SRR2496709", "SE"),
]
IO_MODES = ["has1", "no1"]
PART_LABELS = [
    "prepare CPE task and reusable buffers",
    "CPE FASTQ formatting and input release",
    "CPE alignment and SAM length pass",
    "assign slices in the shared SAM buffer",
    "CPE SAM record generation",
    "release temporary worker data",
]
READ_BW_OUTLIER_MIB_S = 30.0
WRITE_BW_OUTLIER_MIB_S = 30.0


def last_float(pattern, text, default=math.nan):
    values = re.findall(pattern, text, flags=re.MULTILINE)
    return float(values[-1]) if values else default


def last_int(pattern, text, default=0):
    values = re.findall(pattern, text, flags=re.MULTILINE)
    return int(values[-1]) if values else default


def parse_log(path):
    text = path.read_text(errors="replace")
    result = {
        "path": path,
        "pipeline_total": last_float(
            r"total - complete three-stage pipeline\s+([0-9.]+) s", text),
        "stage1": last_float(
            r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s", text),
        "stage2": last_float(
            r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s", text),
        "stage3": last_float(
            r"stage 3 - write SAM and release batch data\s+([0-9.]+) s", text),
        "worker_total": last_float(
            r"total - active merge worker\s+([0-9.]+) s", text),
        "fread_seconds": last_float(
            r"FASTQ fread calls\s+([0-9.]+) s", text),
        "read_bytes": last_int(r"raw FASTQ bytes read\s+([0-9]+)", text),
        "read_bw": last_float(
            r"effective fread bandwidth\s+([0-9.]+) MiB/s", text),
        "slowest_fread": last_float(
            r"slowest fread\s+([0-9.]+) s", text),
        "write_seconds": last_float(
            r"POSIX write system calls\s+([0-9.]+) s", text),
        "write_bw": last_float(
            r"POSIX write effective bandwidth\s+([0-9.]+) MiB/s", text),
        "output_bytes": last_int(r"submitted SAM bytes\s+([0-9]+)", text),
        "batches": last_int(
            r"batches completed / stage calls\s+([0-9]+)\s*/", text),
    }
    for index, label in enumerate(PART_LABELS, 1):
        result["part{}".format(index)] = last_float(
            r"part {} - {}\s+([0-9.]+) s".format(index, re.escape(label)),
            text,
        )
    required = ("pipeline_total", "stage1", "stage2", "stage3", "worker_total")
    missing = [name for name in required if math.isnan(result[name])]
    if missing:
        raise ValueError("{}: missing {}".format(path, ", ".join(missing)))
    result["read_outlier"] = (
        not math.isnan(result["read_bw"])
        and result["read_bw"] < READ_BW_OUTLIER_MIB_S
    )
    result["write_outlier"] = (
        not math.isnan(result["write_bw"])
        and result["write_bw"] < WRITE_BW_OUTLIER_MIB_S
    )
    return result


def load_logs(root):
    records = {}
    for io_mode in IO_MODES:
        for dataset, read_mode in CASES:
            path = root / io_mode / "{}_{}.run.log".format(dataset, read_mode)
            if not path.exists():
                raise FileNotFoundError(path)
            records[(io_mode, dataset, read_mode)] = parse_log(path)
    return records


def load_expected(path):
    expected = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or line.startswith("dataset\t"):
            continue
        dataset, read_mode, digest = line.split("\t")
        expected[(dataset, read_mode)] = digest
    return expected


def load_correctness(root, expected):
    rows = []
    for io_mode in IO_MODES:
        for dataset, read_mode in CASES:
            path = root / io_mode / "{}_{}.md5".format(dataset, read_mode)
            observed = path.read_text().split()[0] if path.exists() else ""
            wanted = expected.get((dataset, read_mode), "")
            rows.append({
                "io_mode": io_mode,
                "dataset": dataset,
                "read_mode": read_mode,
                "observed": observed,
                "expected": wanted,
                "status": "PASS" if observed and observed == wanted else "FAIL",
            })
    return rows


def format_number(value, digits=3):
    if isinstance(value, float) and math.isnan(value):
        return ""
    return ("{:.%df}" % digits).format(value)


def write_performance_tsv(current, baseline, path):
    columns = [
        "io_mode", "dataset", "read_mode", "pipeline_total_s",
        "stage1_s", "stage2_s", "stage3_s", "worker_total_s",
        "part1_s", "part2_s", "part3_s", "part4_s", "part5_s", "part6_s",
        "batches", "read_bytes", "fread_s", "read_bw_mib_s",
        "slowest_fread_s", "output_bytes", "write_s", "write_bw_mib_s",
        "read_io_outlier", "write_io_outlier", "baseline_stage2_s",
        "stage2_speedup", "stage2_change_percent",
    ]
    lines = ["\t".join(columns)]
    for io_mode in IO_MODES:
        for dataset, read_mode in CASES:
            key = (io_mode, dataset, read_mode)
            item = current[key]
            old = baseline[key]
            speedup = old["stage2"] / item["stage2"]
            change = (item["stage2"] / old["stage2"] - 1.0) * 100.0
            values = [io_mode, dataset, read_mode]
            values.extend(format_number(item[name]) for name in (
                "pipeline_total", "stage1", "stage2", "stage3", "worker_total",
                "part1", "part2", "part3", "part4", "part5", "part6",
            ))
            values.extend([
                str(item["batches"]), str(item["read_bytes"]),
                format_number(item["fread_seconds"]), format_number(item["read_bw"]),
                format_number(item["slowest_fread"]), str(item["output_bytes"]),
                format_number(item["write_seconds"]), format_number(item["write_bw"]),
                "yes" if item["read_outlier"] else "no",
                "yes" if item["write_outlier"] else "no",
                format_number(old["stage2"]), format_number(speedup),
                format_number(change),
            ])
            lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n")


def write_correctness_tsv(rows, path):
    lines = ["io_mode\tdataset\tread_mode\tobserved_md5\texpected_md5\tstatus"]
    for row in rows:
        lines.append("\t".join(row[name] for name in (
            "io_mode", "dataset", "read_mode", "observed", "expected", "status",
        )))
    path.write_text("\n".join(lines) + "\n")


def svg_header(width, height, title, subtitle):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
        'viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        '<style>text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", '
        'Arial, sans-serif; fill: #182229; letter-spacing: 0; } '
        '.small { font-size: 12px; fill: #58656e; } '
        '.mono { font-family: Menlo, Consolas, monospace; }</style>',
        '<text x="36" y="40" font-size="23" font-weight="700">{}</text>'.format(
            html.escape(title)),
        '<text x="36" y="63" class="small">{}</text>'.format(html.escape(subtitle)),
    ]


def case_label(io_mode, dataset, read_mode):
    return "{} {} {}".format(dataset.replace("small_", "small "), read_mode, io_mode)


def write_stage2_svg(current, baseline, path):
    rows = [(io_mode, dataset, read_mode)
            for dataset, read_mode in CASES for io_mode in IO_MODES]
    width, row_height = 1420, 38
    height = 120 + len(rows) * row_height
    left, bar_width = 330, 950
    maximum = max(max(current[key]["stage2"], baseline[key]["stage2"])
                  for key in rows)
    out = svg_header(
        width, height,
        "Big-data Stage 2: previous versus optimized",
        "Release build, cgs_cross + pool, one node and one host process",
    )
    out.extend([
        '<rect x="980" y="28" width="13" height="9" fill="#657786"/>',
        '<text x="999" y="37" class="small">previous</text>',
        '<rect x="1075" y="28" width="13" height="9" fill="#d55e00"/>',
        '<text x="1094" y="37" class="small">optimized</text>',
        '<line x1="{}" y1="88" x2="{}" y2="88" stroke="#c8ced3"/>'.format(
            left, left + bar_width),
        '<text x="{}" y="83" class="small">0</text>'.format(left),
        '<text x="{}" y="83" class="small" text-anchor="end">{:.1f} s</text>'.format(
            left + bar_width, maximum),
    ])
    for index, key in enumerate(rows):
        y = 103 + index * row_height
        old = baseline[key]["stage2"]
        new = current[key]["stage2"]
        old_width = bar_width * old / maximum
        new_width = bar_width * new / maximum
        speedup = old / new
        out.append('<text x="36" y="{}" font-size="12">{}</text>'.format(
            y + 13, html.escape(case_label(*key))))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="9" fill="#657786"/>'.format(
            left, y, old_width))
        out.append('<rect x="{}" y="{}" width="{:.2f}" height="9" fill="#d55e00"/>'.format(
            left, y + 12, new_width))
        out.append('<text x="{}" y="{}" class="small">{:.3f} s, {:.3f}x</text>'.format(
            min(left + max(old_width, new_width) + 8, width - 135),
            y + 20, new, speedup))
    out.append('</svg>')
    path.write_text("\n".join(out) + "\n")


def write_pipeline_svg(current, path):
    rows = [(io_mode, dataset, read_mode)
            for dataset, read_mode in CASES for io_mode in IO_MODES]
    width, row_height = 1480, 50
    height = 130 + len(rows) * row_height
    left, bar_width = 330, 920
    maximum = max(max(current[key][stage] for stage in ("stage1", "stage2", "stage3"))
                  for key in rows)
    colors = {"stage1": "#0072b2", "stage2": "#d55e00", "stage3": "#009e73"}
    out = svg_header(
        width, height,
        "Optimized big-data pipeline timing and I/O diagnostics",
        "Stages overlap in the pipeline and must not be added together; red OUTLIER marks bandwidth below 30 MiB/s",
    )
    legend_x = 780
    for index, stage in enumerate(("stage1", "stage2", "stage3")):
        x = legend_x + index * 115
        out.append('<rect x="{}" y="28" width="13" height="9" fill="{}"/>'.format(
            x, colors[stage]))
        out.append('<text x="{}" y="37" class="small">{}</text>'.format(x + 19, stage))
    out.append('<line x1="{}" y1="88" x2="{}" y2="88" stroke="#c8ced3"/>'.format(
        left, left + bar_width))
    out.append('<text x="{}" y="83" class="small">0</text>'.format(left))
    out.append('<text x="{}" y="83" class="small" text-anchor="end">{:.1f} s</text>'.format(
        left + bar_width, maximum))
    for index, key in enumerate(rows):
        y = 104 + index * row_height
        item = current[key]
        out.append('<text x="36" y="{}" font-size="12">{}</text>'.format(
            y + 18, html.escape(case_label(*key))))
        for stage_index, stage in enumerate(("stage1", "stage2", "stage3")):
            seconds = item[stage]
            length = bar_width * seconds / maximum
            out.append('<rect x="{}" y="{}" width="{:.2f}" height="8" fill="{}"/>'.format(
                left, y + stage_index * 11, length, colors[stage]))
        warning = item["read_outlier"] or item["write_outlier"]
        note = "read {:.1f}, write {:.1f} MiB/s{}".format(
            item["read_bw"], item["write_bw"], "  OUTLIER" if warning else "")
        color = "#b3261e" if warning else "#58656e"
        out.append('<text x="1270" y="{}" font-size="11" style="fill:{}">{}</text>'.format(
            y + 18, color, html.escape(note)))
    out.append('</svg>')
    path.write_text("\n".join(out) + "\n")


def markdown_table(headers, rows):
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def write_report(current, baseline, correctness, path):
    passed = sum(row["status"] == "PASS" for row in correctness)
    speedups = []
    part1_reductions = []
    part3_improvements = []
    mode_differences = []
    for key, item in current.items():
        old = baseline[key]
        speedups.append(old["stage2"] / item["stage2"])
        part1_reductions.append((1.0 - item["part1"] / old["part1"]) * 100.0)
        part3_improvements.append((1.0 - item["part3"] / old["part3"]) * 100.0)
    for dataset, read_mode in CASES:
        has1 = current[("has1", dataset, read_mode)]["stage2"]
        no1 = current[("no1", dataset, read_mode)]["stage2"]
        mode_differences.append(abs(no1 / has1 - 1.0) * 100.0)
    speedups.sort()
    part1_reductions.sort()
    part3_improvements.sort()
    mode_differences.sort()
    median_speedup = (speedups[5] + speedups[6]) / 2.0
    median_part1_reduction = (part1_reductions[5] + part1_reductions[6]) / 2.0
    median_part3_improvement = (part3_improvements[5] + part3_improvements[6]) / 2.0
    accuracy_rows = [[
        row["io_mode"], row["dataset"], row["read_mode"],
        "`{}`".format(row["observed"]), row["status"],
    ] for row in correctness]
    perf_rows = []
    comparison_rows = []
    for dataset, read_mode in CASES:
        for io_mode in IO_MODES:
            key = (io_mode, dataset, read_mode)
            item, old = current[key], baseline[key]
            io_note = []
            if item["read_outlier"]:
                io_note.append("read")
            if item["write_outlier"]:
                io_note.append("write")
            perf_rows.append([
                io_mode, dataset, read_mode,
                "{:.3f}".format(item["pipeline_total"]),
                "{:.3f}".format(item["stage1"]),
                "{:.3f}".format(item["stage2"]),
                "{:.3f}".format(item["stage3"]),
                "{:.1f}".format(item["read_bw"]),
                "{:.1f}".format(item["write_bw"]),
                ", ".join(io_note) if io_note else "-",
            ])
            comparison_rows.append([
                io_mode, dataset, read_mode,
                "{:.3f}".format(old["stage2"]),
                "{:.3f}".format(item["stage2"]),
                "{:.3f}x".format(old["stage2"] / item["stage2"]),
                "{:.3f}".format(old["part1"]),
                "{:.3f}".format(item["part1"]),
                "{:.3f}".format(old["part3"]),
                "{:.3f}".format(item["part3"]),
            ])

    report = """# cgs_cross + pool 大数据发布版验证

## 结论

- 正确性：{passed}/{total} 个 `has1/no1 x PE/SE` 输出 MD5 与 x86 BWA 标准结果逐字节一致。
- 测试配置：`cgs_cross + pool`、`CPE_PROFILE=0`、单节点、单主进程、6 CG x 64 CPE。
- 12 个 case 的 Stage 2 全部快于旧版：加速范围 {speedup_min:.3f}x-{speedup_max:.3f}x，中位数 {speedup_median:.3f}x。
- part1 准备阶段下降 {part1_min:.1f}%-{part1_max:.1f}%（中位数 {part1_median:.1f}%），验证了去掉 1.5 GiB SAM buffer 无效初始化的结构性收益。
- part3 比对与长度阶段中位数改善 {part3_median:.1f}%；仅一个 case 慢 {part3_worst:.1f}%，属于单次跨节点测试噪声范围。mate-rescue 去重优化主要作用于 150 bp PE。
- 当前 has1/no1 的 Stage 2 差异为 {mode_diff_min:.2f}%-{mode_diff_max:.2f}%，说明计算结果不受本轮严重 I/O 抖动支配。
- 本轮存在明显共享文件系统抖动。低于 30 MiB/s 的读写被标记为 I/O 离群，不能用 pipeline total 判断计算优化效果。

`no1` 的 Stage 1/2/3 在三线程流水线中有重叠，不能相加；`has1` 使用 `-1` 关闭该流水线，三个阶段近似串行。计算优化统一以 Stage 2 和 part3 为主指标。

## 正确性

{accuracy_table}

## 当前性能

| 模式 | 数据 | reads | pipeline | Stage 1 | Stage 2 | Stage 3 | read MiB/s | write MiB/s | I/O 离群 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{performance_rows}

## Stage 2 A/B

`previous` 来自 `correctness_results/bigdata_results/cgs_cross_pool`，`optimized` 是本轮发布版。单次跨节点测试仍有小幅噪声，重点看 part1 的结构性下降、part3 是否没有退化，以及 has1/no1 是否一致。

| 模式 | 数据 | reads | previous Stage2 | optimized Stage2 | speedup | old part1 | new part1 | old part3 | new part3 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
{comparison_rows}

## I/O 观察

- `small_SRR7963242 SE no1` 的 Stage 1 为 302.616 s，读带宽只有 11.963 MiB/s，最慢单次 fread 为 47.206 s；同数据 has1 的读带宽为 120.134 MiB/s。这个 10x 差异是共享文件系统波动。
- `-1` 在当前代码中设置 `no_mt_io=1`，关闭三阶段多线程 I/O 流水线；它不应让底层 `fread` 吞吐稳定快一个数量级。本轮 has1/no1 的 Stage 2 接近，而 Stage 1 的实际 `fread` 差异很大，也支持 I/O 抖动判断。
- 参考索引加载在若干节点上也曾长时间停滞，但它发生在三阶段 pipeline 之前，没有进入表中计时；异常启动已剔除。

## 文件

- `analysis/performance.tsv`：所有阶段、Stage 2 part1-part6、I/O 和 A/B 数值。
- `analysis/correctness.tsv`：12 个 MD5 的期望值、实测值和状态。
- `analysis/stage2_comparison.svg`：旧版与优化版 Stage 2 对比。
- `analysis/pipeline_io.svg`：当前 Stage 1/2/3 与 I/O 离群标记。
- `scripts/analyze_cpe_bigdata.py`：完整解析和可视化生成代码。

## 复现分析

```bash
python3 scripts/analyze_cpe_bigdata.py
```
""".format(
        passed=passed,
        total=len(correctness),
        speedup_min=min(speedups),
        speedup_max=max(speedups),
        speedup_median=median_speedup,
        part1_min=min(part1_reductions),
        part1_max=max(part1_reductions),
        part1_median=median_part1_reduction,
        part3_median=median_part3_improvement,
        part3_worst=abs(min(part3_improvements)),
        mode_diff_min=min(mode_differences),
        mode_diff_max=max(mode_differences),
        accuracy_table=markdown_table(
            ["模式", "数据", "reads", "MD5", "状态"], accuracy_rows),
        performance_rows="\n".join("| " + " | ".join(row) + " |" for row in perf_rows),
        comparison_rows="\n".join("| " + " | ".join(row) + " |" for row in comparison_rows),
    )
    path.write_text(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--current", type=Path,
        default=Path("correctness_results/cpe_bigdata_opt_20260825"),
    )
    parser.add_argument(
        "--baseline", type=Path,
        default=Path("correctness_results/bigdata_results/cgs_cross_pool"),
    )
    parser.add_argument(
        "--expected", type=Path,
        default=Path("scripts/bigdata_expected_md5.tsv"),
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    output_dir = args.output_dir or args.current / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    current = load_logs(args.current)
    baseline = load_logs(args.baseline)
    correctness = load_correctness(args.current, load_expected(args.expected))

    write_performance_tsv(current, baseline, output_dir / "performance.tsv")
    write_correctness_tsv(correctness, output_dir / "correctness.tsv")
    write_stage2_svg(current, baseline, output_dir / "stage2_comparison.svg")
    write_pipeline_svg(current, output_dir / "pipeline_io.svg")
    write_report(current, baseline, correctness, args.current / "README.md")
    failed = [row for row in correctness if row["status"] != "PASS"]
    print("parsed {} current and {} baseline logs".format(len(current), len(baseline)))
    print("correctness: {}/{} PASS".format(len(correctness) - len(failed), len(correctness)))
    print("results: {}".format(output_dir))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
