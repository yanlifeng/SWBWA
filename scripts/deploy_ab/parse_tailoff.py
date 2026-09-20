#!/usr/bin/env python3
"""Parser for the tail-off experiment: parse small/big4 x A/B x rep discard
logs.

Usage: parse_tailoff.py <log_dir>
Output: writes <log_dir>/tailoff_summary.tsv and prints the A/B comparison and
hash consistency to stdout.
"""
import re
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_discard_hash import parse_log as validate_hash

LOG_RE = re.compile(r"^(small|big4)_([AB])_rep(\d+)\.log$")

PAT_CG = re.compile(r"\[MPI input\].*?chunks=(\d+).*?cg_count=(\d+)")
PAT_WORLD = re.compile(r"MPI rank:\s+(\d+)\s*/\s*(\d+)")
PAT_TOTAL = re.compile(r"total - complete three-stage pipeline\s+([0-9.]+) s")
PAT_STAGE1 = re.compile(r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s")
PAT_STAGE2 = re.compile(r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s")
PAT_STAGE3 = re.compile(r"stage 3 - write SAM and release batch data\s+([0-9.]+) s")
PAT_PART3 = re.compile(r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s")
PAT_CLAIMED = re.compile(r"claimed chunks\s+(\d+)")
PAT_RMA = re.compile(r"RMA ticket path total\s+([0-9.]+) s")
PAT_ACCUM2 = re.compile(r"accumulated stage 2\s+([0-9.]+) s")
PAT_BARRIER = re.compile(r"MPI_Barrier \(wait for other ranks\)\s+([0-9.]+) s")
PAT_IPROBE_CALLS = re.compile(r"MPI_Iprobe calls\s+(\d+)")
PAT_IPROBE_TOTAL = re.compile(r"MPI_Iprobe total\s+([0-9.]+) s")
PAT_BW = re.compile(r"effective fread bandwidth\s+([\d.]+) MiB/s")
PAT_FREAD_SLOW = re.compile(r"FASTQ fread calls\s+[\d.]+\s+s\s+\(\d+ calls, max ([\d.]+) s\)")
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
    m = LOG_RE.match(path.name)
    dataset, cfg, rep = m.group(1), m.group(2), int(m.group(3))

    cg_m = PAT_CG.search(text)
    chunks, cg_count = (int(cg_m.group(1)), int(cg_m.group(2))) if cg_m else (None, None)

    world = 0
    for w in PAT_WORLD.finditer(text):
        world = max(world, int(w.group(2)))

    def grab(pat):
        return [float(x) for x in pat.findall(text)]

    hash_sum = 0
    hash_xor = 0
    hash_ranks = set()
    hash_bytes = 0
    for hm in PAT_HASH.finditer(text):
        hash_ranks.add(int(hm.group(1)))
        hash_sum = (hash_sum + int(hm.group(5), 16)) & MASK64
        hash_xor ^= int(hm.group(6), 16)
        hash_bytes += int(hm.group(4))

    try:
        validate_hash(path)
        hash_valid = True
    except (OSError, ValueError):
        hash_valid = False

    return {
        "log": path.name,
        "dataset": dataset,
        "config": cfg,
        "rep": rep,
        "cg_count": cg_count,
        "mpi_chunks": chunks,
        "world_size": world,
        "hash_ranks": len(hash_ranks),
        "hash_valid": hash_valid,
        "hash_bytes": hash_bytes,
        "hash_sum": f"{hash_sum:016x}",
        "hash_xor": f"{hash_xor:016x}",
        "total": summarize(grab(PAT_TOTAL)),
        "stage1": summarize(grab(PAT_STAGE1)),
        "stage2": summarize(grab(PAT_STAGE2)),
        "stage3": summarize(grab(PAT_STAGE3)),
        "part3": summarize(grab(PAT_PART3)),
        "claimed": summarize(grab(PAT_CLAIMED)),
        "rma": summarize(grab(PAT_RMA)),
        "accum_stage2": summarize(grab(PAT_ACCUM2)),
        "barrier": summarize(grab(PAT_BARRIER)),
        "iprobe_calls": summarize(grab(PAT_IPROBE_CALLS)),
        "iprobe_total": summarize(grab(PAT_IPROBE_TOTAL)),
        "read_bw": summarize(grab(PAT_BW)),
        "fread_slowest": summarize(grab(PAT_FREAD_SLOW)),
    }


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    logs = sorted(p for p in root.iterdir() if LOG_RE.match(p.name))
    if not logs:
        print(f"no matching logs in {root}", file=sys.stderr)
        return 1

    rows = [parse_log(p) for p in logs]

    cols = ["log", "dataset", "config", "rep", "cg_count", "mpi_chunks",
            "world_size", "hash_ranks", "hash_valid", "hash_bytes", "hash_sum", "hash_xor",
            "total", "stage1", "stage2", "stage3", "part3", "claimed",
            "rma", "accum_stage2", "barrier", "iprobe_calls", "iprobe_total",
            "read_bw", "fread_slowest"]

    def cell(r, c):
        v = r[c]
        return fmt(v) if isinstance(v, tuple) else ("" if v is None else str(v))

    out = root / "tailoff_summary.tsv"
    with out.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(cell(r, c) for c in cols) + "\n")

    print("=" * 100)
    print("逐次运行")
    print("=" * 100)
    for r in rows:
        print(f"[{r['dataset']} {r['config']} rep{r['rep']}] chunks={r['mpi_chunks']} "
              f"world={r['world_size']} hash_ranks={r['hash_ranks']}")
        print(f"  total={fmt(r['total'])} stage1={fmt(r['stage1'])} "
              f"stage2={fmt(r['stage2'])} stage3={fmt(r['stage3'])}")
        print(f"  part3={fmt(r['part3'])} rma={fmt(r['rma'])} "
              f"barrier={fmt(r['barrier'])} read_bw={fmt(r['read_bw'])}")

    # A/B comparison (per dataset x metric, median of 3 reps).
    print()
    print("=" * 100)
    print("A/B 对照（median of 3 reps per dataset）")
    print("=" * 100)
    metrics = ["total", "stage1", "stage2", "stage3", "part3",
               "rma", "accum_stage2", "barrier", "read_bw"]
    for ds in ("small", "big4"):
        print(f"\n--- dataset {ds} ---")
        print(f"{'metric':<14}{'A':>22}{'B':>22}{'B/A':>10}")
        for met in metrics:
            a = [r[met][3] for r in rows if r["dataset"] == ds and r["config"] == "A" and r[met][0] > 0]
            b = [r[met][3] for r in rows if r["dataset"] == ds and r["config"] == "B" and r[met][0] > 0]
            am = median(a) if a else float("nan")
            bm = median(b) if b else float("nan")
            ratio = (bm / am) if am and am == am and bm == bm and am > 0 else float("nan")
            print(f"{met:<14}{am:>22.3f}{bm:>22.3f}{ratio:>10.3f}")

    # hash consistency.
    print()
    print("跨配置指纹一致性（sum/xor 均顺序无关，A 与 B 应相等）")
    for ds in ("small", "big4"):
        s = {r["hash_sum"] for r in rows if r["dataset"] == ds}
        x = {r["hash_xor"] for r in rows if r["dataset"] == ds}
        valid = all(r["hash_valid"] for r in rows if r["dataset"] == ds)
        tag = "一致" if (valid and len(s) == 1 and len(x) == 1) else "不一致或未验证!"
        print(f"  {ds}: distinct sum={len(s)} xor={len(x)} -> {tag}  {sorted(s)}")

    print(f"\n写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
