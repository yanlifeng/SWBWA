#!/usr/bin/env python3
"""用 max（最慢进程/makespan）重新统计 tailoff A/B 实验。

MPI 程序的墙钟时间由最慢进程决定，必须用 max 而非 median。
stage1/stage2/stage3/part3/barrier 都是每进程累加值，max 即关键路径。
"""
import re
from pathlib import Path
from statistics import median, mean

LOG_RE = re.compile(r"^(small|big4)_([AB])_rep(\d+)\.log$")
PAT_TOTAL = re.compile(r"total - complete three-stage pipeline\s+([0-9.]+) s")
PAT_STAGE1 = re.compile(r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s")
PAT_STAGE2 = re.compile(r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s")
PAT_STAGE3 = re.compile(r"stage 3 - write SAM and release batch data\s+([0-9.]+) s")
PAT_PART3 = re.compile(r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s")
PAT_BARRIER = re.compile(r"MPI_Barrier \(wait for other ranks\)\s+([0-9.]+) s")
PAT_RMA = re.compile(r"RMA ticket path total\s+([0-9.]+) s")
PAT_ACCUM2 = re.compile(r"accumulated stage 2\s+([0-9.]+) s")
PAT_CLAIMED = re.compile(r"claimed chunks\s+(\d+)")


def grab(pat, text):
    return [float(x) for x in pat.findall(text)]


def stat(vals):
    if not vals:
        return None
    return {
        "n": len(vals),
        "max": max(vals),
        "min": min(vals),
        "mean": mean(vals),
        "median": median(vals),
    }


rows = []
for p in sorted(Path(".").glob("*.log")):
    m = LOG_RE.match(p.name)
    if not m:
        continue
    ds, cfg, rep = m.group(1), m.group(2), int(m.group(3))
    t = p.read_text(errors="replace")
    rows.append({
        "ds": ds, "cfg": cfg, "rep": rep,
        "total": stat(grab(PAT_TOTAL, t)),
        "stage1": stat(grab(PAT_STAGE1, t)),
        "stage2": stat(grab(PAT_STAGE2, t)),
        "stage3": stat(grab(PAT_STAGE3, t)),
        "part3": stat(grab(PAT_PART3, t)),
        "barrier": stat(grab(PAT_BARRIER, t)),
        "rma": stat(grab(PAT_RMA, t)),
        "accum2": stat(grab(PAT_ACCUM2, t)),
        "claimed": stat(grab(PAT_CLAIMED, t)),
    })

# 每个指标取 3 次 rep 的 max 中位数（makespan 的稳定估计）
def rep_median_max(rows, ds, cfg, key):
    vals = [r[key]["max"] for r in rows if r["ds"] == ds and r["cfg"] == cfg and r[key]]
    return median(vals) if vals else float("nan")


print("=" * 100)
print("makespan 统计（每个指标取 3 次 rep 的 max 中位数，单位：秒）")
print("=" * 100)

metrics = ["total", "stage1", "stage2", "stage3", "part3", "barrier", "rma"]
for ds in ("small", "big4"):
    print(f"\n--- {ds} ---")
    print(f"{'metric':<10}{'A max':>12}{'B max':>12}{'B/A':>10}")
    for key in metrics:
        a = rep_median_max(rows, ds, "A", key)
        b = rep_median_max(rows, ds, "B", key)
        ratio = b / a if a and a == a and b == b and a > 0 else float("nan")
        print(f"{key:<10}{a:>12.3f}{b:>12.3f}{ratio:>10.3f}")

# 逐次明细
print("\n" + "=" * 100)
print("逐次明细（每行是该 rep 的 max / min / median，n 为采集到的进程数）")
print("=" * 100)
for r in rows:
    def f(key):
        s = r[key]
        return f"max={s['max']:.3f} min={s['min']:.3f} med={s['median']:.3f} (n={s['n']})" if s else "n/a"
    print(f"[{r['ds']} {r['cfg']} rep{r['rep']}]")
    print(f"  total  : {f('total')}")
    print(f"  stage1 : {f('stage1')}")
    print(f"  stage2 : {f('stage2')}")
    print(f"  stage3 : {f('stage3')}")
    print(f"  part3  : {f('part3')}")
    print(f"  barrier: {f('barrier')}")
    print(f"  claimed: {f('claimed')}")
