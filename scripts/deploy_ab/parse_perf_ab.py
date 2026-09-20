#!/usr/bin/env python3
"""Fault-tolerant parser for the selling-point-B performance comparison
(A = bigshare 4 procs / B = single-CG 24 procs, discard mode).

Why not reuse scripts/analyze_mpi_discard_profile.py directly:
  that script requires exactly world_size samples per metric, otherwise it
  raises. Cross-node stderr interleaving truncates/drops some ranks' detail
  lines (worse for config B with 24 ranks), making it hard-fail. This script
  instead samples whatever is available and reports count / min / max /
  median. Because samples may be missing, conclusions rely only on the
  min-max range and the median, never on per-rank labels.

Usage: parse_perf_ab.py <perf_log_dir>
Output: writes perf_summary.tsv into that directory and prints the A/B
comparison table to stdout.
"""
import re
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_discard_hash import parse_log as validate_hash

LOG_RE = re.compile(r"^([AB])_.*_(PE|SE)\.run(\d+)\.log$")

PAT_CG = re.compile(r"\[MPI input\].*?chunks=(\d+).*?cg_count=(\d+)")
PAT_WORLD = re.compile(r"MPI rank:\s+(\d+)\s*/\s*(\d+)")
PAT_BWT = re.compile(r"malloc bwt->bwt")

PAT_TOTAL = re.compile(r"total - complete three-stage pipeline\s+([0-9.]+) s")
PAT_STAGE1 = re.compile(r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s")
PAT_STAGE2 = re.compile(r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s")
PAT_STAGE3 = re.compile(r"stage 3 - write SAM and release batch data\s+([0-9.]+) s")

PAT_PART = {
    "part1": re.compile(r"part 1 - prepare CPE task and reusable buffers\s+([0-9.]+) s"),
    "part2": re.compile(r"part 2 - CPE FASTQ formatting and input release\s+([0-9.]+) s"),
    "part3": re.compile(r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s"),
    "part4": re.compile(r"part 4 - assign slices in the shared SAM buffer\s+([0-9.]+) s"),
    "part5": re.compile(r"part 5 - CPE SAM record generation\s+([0-9.]+) s"),
    "part6": re.compile(r"part 6 - release temporary worker data\s+([0-9.]+) s"),
}
PAT_CLAIMED = re.compile(r"claimed chunks\s+(\d+)")
PAT_RMA = re.compile(r"RMA ticket path total\s+([0-9.]+) s")
PAT_ACCUM2 = re.compile(r"accumulated stage 2\s+([0-9.]+) s")
PAT_IPROBE_CALLS = re.compile(r"MPI_Iprobe calls\s+(\d+)")
PAT_IPROBE_TOTAL = re.compile(r"MPI_Iprobe total\s+([0-9.]+) s")
PAT_HASH = re.compile(
    r"\[SWBWA output hash rank\s+(\d+)/(\d+)\]\s+calls=(\d+)\s+bytes=(\d+)\s+"
    r"sum=0x([0-9a-fA-F]{16})\s+xor=0x([0-9a-fA-F]{16})(?:\s+enabled=([01]))?"
)
MASK64 = (1 << 64) - 1


def summarize(values):
    if not values:
        return (0, None, None, None)
    return (len(values), min(values), max(values), median(values))


def fmt(t):
    n, lo, hi, med = t
    if n == 0:
        return "n=0"
    return f"n={n} min={lo:.3f} max={hi:.3f} med={med:.3f}"


def parse_log(path):
    text = path.read_text(errors="replace")
    name = LOG_RE.match(path.name)
    cfg, read_mode, rep = name.groups()

    cg_m = PAT_CG.search(text)
    chunks, cg_count = (int(cg_m.group(1)), int(cg_m.group(2))) if cg_m else (None, None)

    world = 0
    for m in PAT_WORLD.finditer(text):
        world = max(world, int(m.group(2)))

    def grab(pat):
        return [float(x) for x in pat.findall(text)]

    hash_sum = 0
    hash_xor = 0
    hash_ranks = set()
    hash_enabled = None
    for m in PAT_HASH.finditer(text):
        rank = int(m.group(1))
        hash_ranks.add(rank)
        hash_sum = (hash_sum + int(m.group(5), 16)) & MASK64
        hash_xor ^= int(m.group(6), 16)
        if m.group(7) is not None:
            hash_enabled = m.group(7)

    try:
        validate_hash(path)
        hash_valid = True
    except (OSError, ValueError):
        hash_valid = False

    return {
        "log": path.name,
        "config": cfg,
        "read_mode": read_mode,
        "rep": int(rep),
        "cg_count": cg_count,
        "mpi_chunks": chunks,
        "world_size": world,
        "world_size_observed": world,
        "hash_ranks": len(hash_ranks),
        "hash_enabled": hash_enabled,
        "hash_valid": hash_valid,
        "hash_sum": f"{hash_sum:016x}",
        "hash_xor": f"{hash_xor:016x}",
        "bwt_copies": len(PAT_BWT.findall(text)),
        "total": summarize(grab(PAT_TOTAL)),
        "stage1": summarize(grab(PAT_STAGE1)),
        "stage2": summarize(grab(PAT_STAGE2)),
        "stage3": summarize(grab(PAT_STAGE3)),
        "part3": summarize(grab(PAT_PART["part3"])),
        "claimed": summarize(grab(PAT_CLAIMED)),
        "rma": summarize(grab(PAT_RMA)),
        "accum_stage2": summarize(grab(PAT_ACCUM2)),
        "iprobe_calls": summarize(grab(PAT_IPROBE_CALLS)),
        "iprobe_total": summarize(grab(PAT_IPROBE_TOTAL)),
    }


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    logs = sorted(p for p in root.iterdir() if LOG_RE.match(p.name))
    if not logs:
        print(f"no matching logs in {root}", file=sys.stderr)
        return 1

    rows = [parse_log(p) for p in logs]

    cols = ["log", "config", "rep", "cg_count", "mpi_chunks", "world_size",
            "hash_ranks", "hash_enabled", "hash_valid", "hash_sum", "hash_xor", "bwt_copies",
            "total", "stage1", "stage2", "stage3", "part3", "claimed",
            "rma", "accum_stage2", "iprobe_calls", "iprobe_total"]

    def cell(r, c):
        v = r[c]
        if isinstance(v, tuple):
            return fmt(v)
        return "" if v is None else str(v)

    out = root / "perf_summary.tsv"
    with out.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(cell(r, c) for c in cols) + "\n")

    # print key metrics per run.
    print("=" * 96)
    print("逐次运行（跨 rank 的 min/max/median）")
    print("=" * 96)
    for r in rows:
        print(f"\n[{r['config']} run{r['rep']:02d}] {r['log']}")
        print(f"  cg_count={r['cg_count']}  mpi_chunks={r['mpi_chunks']}  "
              f"world_size={r['world_size']}  bwt_copies={r['bwt_copies']}  "
              f"hash_ranks={r['hash_ranks']} enabled={r['hash_enabled']}")
        print(f"  total    : {fmt(r['total'])}")
        print(f"  stage1   : {fmt(r['stage1'])}")
        print(f"  stage2   : {fmt(r['stage2'])}")
        print(f"  stage3   : {fmt(r['stage3'])}")
        print(f"  part3    : {fmt(r['part3'])}")
        print(f"  rma      : {fmt(r['rma'])}")
        print(f"  iprobe   : {fmt(r['iprobe_calls'])} ")
        print(f"  hash_sum={r['hash_sum']} hash_xor={r['hash_xor']}")

    # A/B comparison.
    print()
    print("=" * 96)
    print("A/B 对照（各指标取该配置 3 次运行的中位数）")
    print("=" * 96)
    metrics = ["total", "stage1", "stage2", "stage3", "part3", "rma", "accum_stage2"]
    print(f"{'metric':<16}{'A (4 procs)':>18}{'B (24 procs)':>18}{'B/A':>10}")
    for m in metrics:
        a = [r[m][3] for r in rows if r["config"] == "A" and r[m][0] > 0]
        b = [r[m][3] for r in rows if r["config"] == "B" and r[m][0] > 0]
        am = median(a) if a else float("nan")
        bm = median(b) if b else float("nan")
        ratio = (bm / am) if am and am == am and bm == bm and am > 0 else float("nan")
        print(f"{m:<16}{am:>18.3f}{bm:>18.3f}{ratio:>10.3f}")

    print()
    print("跨配置指纹一致性（sum 与 xor 均与记录顺序无关，应相等）")
    for cfg in ("A", "B"):
        s = {r["hash_sum"] for r in rows if r["config"] == cfg}
        x = {r["hash_xor"] for r in rows if r["config"] == cfg}
        print(f"  {cfg}: distinct hash_sum={len(s)} distinct hash_xor={len(x)}  {sorted(s)}")
    allsum = {r["hash_sum"] for r in rows}
    allxor = {r["hash_xor"] for r in rows}
    print(f"  全部运行: distinct sum={len(allsum)}  distinct xor={len(allxor)}")
    valid = all(r["hash_valid"] for r in rows)
    print(f"  完整全量指纹可用 -> {valid}")
    print(f"  指纹一致 -> {'YES' if valid and len(allsum) == 1 and len(allxor) == 1 else 'NO / UNVERIFIED'}")

    print(f"\n写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
