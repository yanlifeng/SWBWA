#!/usr/bin/env python3
"""Parse the A/B + LWPF-profile run logs produced by ab_cond.sh / ab_driver.sh.

Layout:
    <logs_root>/<cfg>_<dataset>_<mode>/{old_np,new_np,new_p,old_p}.run.log

Emits four TSVs (or to stdout when out_dir is omitted):
  * ab_stage2.tsv              - paired old_np vs new_np stage1/2/3/total
  * kernel_cycles.tsv          - per-kernel LWPF CYCLE average for new_p / old_p
  * kernel_diff.tsv            - the 9 regions both versions instrument, delta + share
  * kernel_new_breakdown.tsv   - all 33 regions of the new build, share of CPE total

Usage: parse_ab.py <logs_root> [out_dir]
"""
import os
import re
import sys

TOTAL_RE = re.compile(r"total - complete three-stage pipeline\s+([0-9.]+)\s*s")
# NB: the logs are CRLF and the timing lines are mid-file, so MULTILINE is
# required for '$' to anchor at a line end rather than end-of-string.
STAGE_RE = re.compile(r"stage ([0-9]) - .*?\s+([0-9.]+)\s*s\s*$", re.MULTILINE)
FREAD_RE = re.compile(r"effective fread bandwidth\s+([0-9.]+)\s*MiB/s")
NUM_RE = re.compile(r"^\s*([0-9.]+)\s*([KMG]?)\.?\s*$")
# CYCLE column of a kernel table row, e.g. "7343M", "  97M", "   0."
CYC_COL_RE = re.compile(r"^\s*[0-9][0-9.]*\s*[KMG]?\.?\s*$")

# Order of LWPF_KERNELS. Verified against the sources:
#   old = fc658b8:slave/swbwa_cpe_profile.h  -> exactly KERNELS[0:9]
#   new = HEAD:src/slave/swbwa_cpe_profile.h -> all 33
# The old 9 are a prefix of the new 33, with identical names and order, so a
# single positional mapping is correct for both variants. Rows 10..33 of the new
# build are sub-regions nested inside rows 1..9 (e.g. MEM_CHAIN ==
# MEM_CHAIN_COLLECT + MEM_CHAIN_BUILD), so rows 1..9 stay comparable across
# versions.
KERNELS = [
    "WORKER_ALIGNMENT", "MEM_CHAIN", "CHAIN_FILTER", "CHAIN_EXTENSION",
    "ALIGNMENT_FINALIZE", "MATE_RESCUE", "PAIRING", "SAM_FORMAT", "SAM_COPY",
    "MEM_CHAIN_COLLECT", "MEM_CHAIN_BUILD", "CHAIN_EXTENSION_DP",
    "MATE_REF_FETCH", "MATE_KSW_ALIGN", "KSW_QUERY_INIT_FORWARD",
    "KSW_DP_FORWARD", "KSW_QUERY_INIT_REVERSE", "KSW_DP_REVERSE", "MATE_DEDUP",
    "DEDUP_SORT_END", "DEDUP_REDUNDANCY", "DEDUP_SORT_SCORE",
    "MEM_COLLECT_FIRST", "MEM_COLLECT_SPLIT", "MEM_COLLECT_LAST",
    "MEM_COLLECT_SORT", "CHAIN_BUILD_REPETITIVE", "CHAIN_BUILD_SA",
    "CHAIN_BUILD_RID", "CHAIN_BUILD_TREE_SEARCH", "CHAIN_BUILD_MERGE",
    "CHAIN_BUILD_INSERT", "CHAIN_BUILD_FINALIZE",
]
# Regions present in the old (fc658b8) build -> the directly comparable set.
OLD_REGION_COUNT = 9
# Only WORKER_ALIGNMENT and SAM_COPY are top-level. SAM_FORMAT is NESTED inside
# WORKER_ALIGNMENT in both modes: slave.c opens WORKER_ALIGNMENT, then
# worker12_pre_fast opens/closes SAM_FORMAT around mem_sam_pe, then
# WORKER_ALIGNMENT closes (see worker12_s_pre_fast / _cross). SAM_COPY is a
# separate task function, so it is the only true sibling. Using
# WORKER_ALIGNMENT + SAM_COPY as the CPE total avoids double counting
# SAM_FORMAT (and MATE_RESCUE/PAIRING, which sit inside SAM_FORMAT in PE mode).
TOP_LEVEL = ("WORKER_ALIGNMENT", "SAM_COPY")
SCALE = {"": 1.0, "K": 1e3, "M": 1e6, "G": 1e9}


def to_num(tok):
    m = NUM_RE.match(tok)
    if not m:
        return float("nan")
    return float(m.group(1)) * SCALE[m.group(2)]


def read_log(path):
    try:
        with open(path, "r", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None

    out = {}
    m = TOTAL_RE.search(text)
    if m:
        out["total"] = float(m.group(1))
        for sm in STAGE_RE.finditer(text):
            out["stage" + sm.group(1)] = float(sm.group(2))
    fm = FREAD_RE.search(text)
    if fm:
        out["fread_bw"] = float(fm.group(1))

    # LWPF kernel table.
    cycles = []
    in_table = False
    for line in text.splitlines():
        if "LWPF kernel summary" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("========"):
            break
        if not line.startswith("|"):
            continue
        fields = line.split("|")
        if len(fields) < 4:
            continue
        name = fields[1].strip().strip(".")
        # Keep only genuine data rows: a (truncated) name plus a numeric CYCLE.
        if not name or name in ("KERNEL", "CYC"):
            continue
        if not CYC_COL_RE.match(fields[2]):
            continue
        cycles.append((name, to_num(fields[2])))
    if cycles:
        out["cycles"] = cycles
    return out if out else None


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: parse_ab.py <logs_root> [out_dir]")
    root = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else None

    stage_rows = ["cfg\tdataset\tmode\told_np_s2\tnew_np_s2\td_s2_pct"
                  "\told_np_total\tnew_np_total\td_total_pct"
                  "\told_np_s1\tnew_np_s1\told_np_s3\tnew_np_s3"
                  "\told_np_bw\tnew_np_bw"]
    kernel_rows = ["cfg\tdataset\tmode\tvariant\tkernel\tcycle_avg"]
    diff_rows = ["cfg\tdataset\tmode\tkernel\told_cyc\tnew_cyc\td_pct"
                 "\tshare_of_new_total_pct"]
    bd_rows = ["cfg\tdataset\tmode\tkernel\tnew_cyc\tshare_of_new_total_pct"]

    for entry in sorted(os.listdir(root)):
        d = os.path.join(root, entry)
        if not os.path.isdir(d):
            continue
        parts = entry.split("_")
        if len(parts) < 3:
            continue
        cfg, ds, mode = parts[0], "_".join(parts[1:-1]), parts[-1]

        logs = {lab: read_log(os.path.join(d, lab + ".run.log"))
                for lab in ("old_np", "new_np", "new_p", "old_p")}
        o, n = logs["old_np"], logs["new_np"]
        if o and n and "stage2" in o and "stage2" in n:
            stage_rows.append(
                "%s\t%s\t%s\t%.2f\t%.2f\t%+.1f\t%.2f\t%.2f\t%+.1f"
                "\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f" % (
                    cfg, ds, mode, o["stage2"], n["stage2"],
                    100.0 * (n["stage2"] - o["stage2"]) / o["stage2"],
                    o["total"], n["total"],
                    100.0 * (n["total"] - o["total"]) / o["total"],
                    o.get("stage1", float("nan")), n.get("stage1", float("nan")),
                    o.get("stage3", float("nan")), n.get("stage3", float("nan")),
                    o.get("fread_bw", float("nan")),
                    n.get("fread_bw", float("nan"))))
        for variant in ("new_p", "old_p"):
            rec = logs[variant]
            if not rec or "cycles" not in rec:
                continue
            for i, (name, cyc) in enumerate(rec["cycles"]):
                label = KERNELS[i] if i < len(KERNELS) else "UNKNOWN_%d" % i
                kernel_rows.append("%s\t%s\t%s\t%s\t%s\t%.0f" % (
                    cfg, ds, mode, variant, label, cyc))

        # Paired per-kernel A/B over the regions both versions instrument.
        np_, op = logs["new_p"], logs["old_p"]
        if np_ and op and "cycles" in np_ and "cycles" in op:
            newmap = {KERNELS[i] if i < len(KERNELS) else "UNKNOWN_%d" % i:
                      c for i, (_, c) in enumerate(np_["cycles"])}
            oldmap = {KERNELS[i] if i < len(KERNELS) else "UNKNOWN_%d" % i:
                      c for i, (_, c) in enumerate(op["cycles"])}
            denom = sum(newmap.get(k, 0.0) for k in TOP_LEVEL) or float("nan")
            for k in KERNELS[:OLD_REGION_COUNT]:
                if k not in newmap or k not in oldmap:
                    continue
                a, b = oldmap[k], newmap[k]
                d = 100.0 * (b - a) / a if a else float("nan")
                share = 100.0 * b / denom if denom else float("nan")
                diff_rows.append("%s\t%s\t%s\t%s\t%.3e\t%.3e\t%+.1f\t%.1f" % (
                    cfg, ds, mode, k, a, b, d, share))

        # Full 33-region breakdown of the new build: where the time goes now.
        if np_ and "cycles" in np_:
            nmap = {KERNELS[i] if i < len(KERNELS) else "U%d" % i: c
                    for i, (_, c) in enumerate(np_["cycles"])}
            tot = sum(nmap.get(k, 0.0) for k in TOP_LEVEL) or float("nan")
            for i, (_, c) in enumerate(np_["cycles"]):
                label = KERNELS[i] if i < len(KERNELS) else "UNKNOWN_%d" % i
                bd_rows.append("%s\t%s\t%s\t%s\t%.3e\t%.1f" % (
                    cfg, ds, mode, label, c,
                    100.0 * c / tot if tot else float("nan")))

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "ab_stage2.tsv"), "w") as fh:
            fh.write("\n".join(stage_rows) + "\n")
        with open(os.path.join(out_dir, "kernel_cycles.tsv"), "w") as fh:
            fh.write("\n".join(kernel_rows) + "\n")
        with open(os.path.join(out_dir, "kernel_diff.tsv"), "w") as fh:
            fh.write("\n".join(diff_rows) + "\n")
        with open(os.path.join(out_dir, "kernel_new_breakdown.tsv"), "w") as fh:
            fh.write("\n".join(bd_rows) + "\n")
        print("wrote %s/ab_stage2.tsv (%d rows)" % (out_dir, len(stage_rows) - 1))
        print("wrote %s/kernel_cycles.tsv (%d rows)" % (out_dir, len(kernel_rows) - 1))
        print("wrote %s/kernel_diff.tsv (%d rows)" % (out_dir, len(diff_rows) - 1))
        print("wrote %s/kernel_new_breakdown.tsv (%d rows)" % (out_dir, len(bd_rows) - 1))
    else:
        print("\n".join(stage_rows))
        print()
        print("\n".join(kernel_rows))


if __name__ == "__main__":
    main()
