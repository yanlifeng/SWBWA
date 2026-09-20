#!/usr/bin/env python3
"""Render the A/B + kernel-profiling results as markdown tables."""
import os
import sys
from collections import defaultdict

OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ab/out"
CFG = sys.argv[2] if len(sys.argv) > 2 else "cross"


def rows(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        out = []
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            if len(vals) == len(head):
                out.append(dict(zip(head, vals)))
    return out


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fmt_cyc(v):
    """Cycles in G (1e9) with 1 decimal, matching the log's own resolution."""
    return "%.1f" % (v / 1e9)


stage = [r for r in rows(os.path.join(OUT, "ab_stage2.tsv")) if r["cfg"] == CFG]
diff = [r for r in rows(os.path.join(OUT, "kernel_diff.tsv")) if r["cfg"] == CFG]
bd = [r for r in rows(os.path.join(OUT, "kernel_new_breakdown.tsv")) if r["cfg"] == CFG]
conds = sorted({(r["dataset"], r["mode"]) for r in stage})

print("## %s: stage 2 wall-clock (has1, no profiler)\n" % CFG)
print("| dataset | mode | old (s) | new (s) | delta |")
print("| --- | --- | ---: | ---: | ---: |")
deltas = []
for ds, mode in conds:
    r = next(x for x in stage if x["dataset"] == ds and x["mode"] == mode)
    d = fnum(r["d_s2_pct"])
    deltas.append(d)
    print("| %s | %s | %.2f | %.2f | %+.1f%% |" %
          (ds, mode, fnum(r["old_np_s2"]), fnum(r["new_np_s2"]), d))
print("| **mean** | | | | **%+.1f%%** |" % (sum(deltas) / len(deltas)))
print()

print("## %s: CPE cycles per region, old vs new (G = 1e9 cycles)\n" % CFG)
print("Share is the new build's share of the CPE total")
print("(WORKER_ALIGNMENT + SAM_FORMAT + SAM_COPY).\n")
kernels = []
for r in diff:
    if r["kernel"] not in kernels:
        kernels.append(r["kernel"])
hdr = "| region | " + " | ".join("%s/%s" % c for c in conds) + " | share(new) |"
print(hdr)
print("| --- | " + " | ".join("---:" for _ in conds) + " | ---: |")
for k in kernels:
    cells = []
    for ds, mode in conds:
        r = next((x for x in diff if x["kernel"] == k and x["dataset"] == ds
                  and x["mode"] == mode), None)
        if not r:
            cells.append("-")
            continue
        a, b = fnum(r["old_cyc"]), fnum(r["new_cyc"])
        cells.append("%s>%s (%+.0f%%)" % (fmt_cyc(a), fmt_cyc(b), fnum(r["d_pct"]))
                     if a else "-")
    shares = [fnum(x["share_of_new_total_pct"]) for x in diff if x["kernel"] == k]
    sh = "%s-%s%%" % ("%.0f" % min(shares), "%.0f" % max(shares)) if shares else "-"
    print("| %s | %s | %s |" % (k, " | ".join(cells), sh))
print()

print("## %s: where the time goes now (new build, all 33 regions)\n" % CFG)
print("| region | " + " | ".join("%s/%s" % c for c in conds) + " |")
print("| --- | " + " | ".join("---:" for _ in conds) + " |")
order = []
for r in bd:
    if r["kernel"] not in order:
        order.append(r["kernel"])
order.sort(key=lambda k: -max([fnum(x["share_of_new_total_pct"])
                               for x in bd if x["kernel"] == k] or [0]))
for k in order:
    cells = []
    for ds, mode in conds:
        r = next((x for x in bd if x["kernel"] == k and x["dataset"] == ds
                  and x["mode"] == mode), None)
        cells.append("%.1f%%" % fnum(r["share_of_new_total_pct"]) if r else "-")
    print("| %s | %s |" % (k, " | ".join(cells)))
