#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


LOG_RE = re.compile(r"(.+)_(PE|SE)\.run(\d+)\.log$")
RANK_RE = re.compile(r"MPI rank:\s+(\d+)\s+/\s+(\d+)")
HASH_RE = re.compile(
    r"\[SWBWA output hash rank\s+(\d+)/(\d+)\] "
    r"calls=(\d+) bytes=(\d+) sum=0x([0-9a-fA-F]{16}) "
    r"xor=0x([0-9a-fA-F]{16})(?: enabled=([01]))?"
)
CHUNK_RE = re.compile(
    r"^\s+(\d+)\s+(\d+)\s+(\d+)\s+"
    r"\[\s*(\d+),\s*(\d+)\)\s+(\d+)\s+(\d+)\s+([0-9.]+)"
)

TIMING_PATTERNS = {
    "stage1_seconds": re.compile(
        r"stage 1 - allocate and read raw FASTQ blocks\s+([0-9.]+) s"
    ),
    "stage2_seconds": re.compile(
        r"stage 2 - align reads and generate SAM records\s+([0-9.]+) s"
    ),
    "cpe_part3_seconds": re.compile(
        r"part 3 - CPE alignment and SAM length pass\s+([0-9.]+) s"
    ),
    "fread_seconds": re.compile(r"FASTQ fread calls\s+([0-9.]+) s"),
    "fread_mib_per_second": re.compile(
        r"effective fread bandwidth\s+([0-9.]+) MiB/s"
    ),
}

SCHEDULER_PATTERNS = {
    "claimed_chunks": re.compile(r"claimed chunks\s+(\d+)"),
    "scheduler_rma_seconds": re.compile(
        r"RMA ticket path total\s+([0-9.]+) s"
    ),
    "scheduler_stage2_seconds": re.compile(
        r"accumulated stage 2\s+([0-9.]+) s"
    ),
}

SLOWEST_CHUNK_RE = re.compile(r"slowest chunk\s+(-?\d+)\s*$")
SLOWEST_CHUNK_SECONDS_RE = re.compile(
    r"slowest chunk stage 2\s+([0-9.]+) s"
)


def parse_log(path):
    name_match = LOG_RE.match(path.name)
    if name_match is None:
        raise ValueError(f"unexpected log name: {path.name}")
    dataset, read_mode, repeat = name_match.groups()
    ranks = {}
    chunks = []
    samples = {
        key: [] for key in (*TIMING_PATTERNS, *SCHEDULER_PATTERNS)
    }
    samples["slowest_chunks"] = []
    pending_slowest_chunk = None
    section = None
    rank = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            hash_match = HASH_RE.search(line)
            if hash_match is not None:
                hash_rank = int(hash_match.group(1))
                hash_world_size = int(hash_match.group(2))
                if (hash_world_size <= 0 or hash_rank < 0 or
                        hash_rank >= hash_world_size):
                    raise ValueError(
                        f"{path}: invalid rank/world size "
                        f"{hash_rank}/{hash_world_size}"
                    )
                entry = ranks.setdefault(hash_rank, {})
                if ("world_size" in entry and
                        entry["world_size"] != hash_world_size):
                    raise ValueError(f"{path}: inconsistent MPI world size")
                entry.update(
                    world_size=hash_world_size,
                    hash_calls=int(hash_match.group(3)),
                    hash_bytes=int(hash_match.group(4)),
                    hash_sum=int(hash_match.group(5), 16),
                    hash_xor=int(hash_match.group(6), 16),
                    hash_enabled=(
                        int(hash_match.group(7))
                        if hash_match.group(7) is not None else 1
                    ),
                )

            if "SWBWA Timing Report" in line:
                section = "timing"
                rank = None
                continue
            if "MPI Dynamic Scheduler Debug" in line:
                section = "scheduler"
                rank = None
                continue

            rank_match = RANK_RE.search(line)
            if rank_match is not None and section is not None:
                rank = int(rank_match.group(1))
                ranks.setdefault(rank, {})["world_size"] = int(
                    rank_match.group(2)
                )
                continue
            if rank is None:
                entry = None
            else:
                entry = ranks.setdefault(rank, {})

            chunk_match = CHUNK_RE.match(line)
            if chunk_match is not None:
                order, chunk_id, queue, start, end, byte_count, records, seconds = (
                    chunk_match.groups()
                )
                chunks.append(
                    {
                        "dataset": dataset,
                        "read_mode": read_mode,
                        "repeat": int(repeat),
                        "rank": (
                            rank
                            if section == "scheduler" and rank is not None
                            else -1
                        ),
                        "order": int(order),
                        "chunk_id": int(chunk_id),
                        "queue": int(queue),
                        "start": int(start),
                        "end": int(end),
                        "bytes": int(byte_count),
                        "records": int(records),
                        "stage2_seconds": float(seconds),
                    }
                )
                continue

            timing_match = None
            for key, pattern in TIMING_PATTERNS.items():
                match = pattern.search(line)
                if match is not None:
                    value = float(match.group(1))
                    samples[key].append(value)
                    timing_match = (key, value)
                    break
            if timing_match is not None:
                if section == "timing" and entry is not None:
                    key, value = timing_match
                    entry[key] = value
                continue

            scheduler_match = None
            for key, pattern in SCHEDULER_PATTERNS.items():
                match = pattern.search(line)
                if match is not None:
                    value = (
                        int(match.group(1))
                        if key == "claimed_chunks"
                        else float(match.group(1))
                    )
                    samples[key].append(value)
                    scheduler_match = (key, value)
                    break
            if scheduler_match is not None:
                if section == "scheduler" and entry is not None:
                    key, value = scheduler_match
                    entry[key] = value
                continue

            slowest_match = SLOWEST_CHUNK_RE.search(line)
            if slowest_match is not None:
                pending_slowest_chunk = int(slowest_match.group(1))
                continue
            slowest_seconds_match = SLOWEST_CHUNK_SECONDS_RE.search(line)
            if slowest_seconds_match is not None:
                samples["slowest_chunks"].append(
                    (
                        pending_slowest_chunk,
                        float(slowest_seconds_match.group(1)),
                    )
                )
                pending_slowest_chunk = None
                continue

    return dataset, read_mode, int(repeat), ranks, chunks, samples


def write_tsv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def range_text(values, digits=3):
    if not values:
        return "-"
    return f"{min(values):.{digits}f}-{max(values):.{digits}f}"


def ratio(values):
    return max(values) / min(values) if values and min(values) > 0.0 else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    args = parser.parse_args()

    rank_rows = []
    chunk_rows = []
    run_rows = []
    rma_upper_bounds = []
    for log in sorted(args.result_root.glob("*.run*.log")):
        dataset, read_mode, repeat, ranks, chunks, samples = parse_log(log)
        if not ranks:
            raise ValueError(f"no rank data in {log}")
        chunk_ids = {}
        for chunk in chunks:
            chunk_ids.setdefault(chunk["rank"], []).append(chunk["chunk_id"])
        for rank, entry in sorted(ranks.items()):
            rank_rows.append(
                {
                    "log": log.name,
                    "dataset": dataset,
                    "read_mode": read_mode,
                    "repeat": repeat,
                    "rank": rank,
                    "stage1_seconds": entry.get("stage1_seconds", ""),
                    "stage2_seconds": entry.get("stage2_seconds", ""),
                    "cpe_part3_seconds": entry.get("cpe_part3_seconds", ""),
                    "fread_seconds": entry.get("fread_seconds", ""),
                    "fread_mib_per_second": entry.get(
                        "fread_mib_per_second", ""
                    ),
                    "claimed_chunks": entry.get("claimed_chunks", ""),
                    "chunk_ids": ",".join(
                        str(item) for item in chunk_ids.get(rank, [])
                    ),
                    "scheduler_rma_seconds": entry.get(
                        "scheduler_rma_seconds", ""
                    ),
                    "scheduler_stage2_seconds": entry.get(
                        "scheduler_stage2_seconds", ""
                    ),
                    "hash_calls": entry.get("hash_calls", ""),
                    "hash_bytes": entry.get("hash_bytes", ""),
                    "hash_sum": (
                        f"{entry['hash_sum']:016x}" if "hash_sum" in entry else ""
                    ),
                    "hash_xor": (
                        f"{entry['hash_xor']:016x}" if "hash_xor" in entry else ""
                    ),
                    "hash_enabled": entry.get("hash_enabled", ""),
                }
            )

        complete = [entry for entry in ranks.values() if "hash_sum" in entry]
        world_sizes = {
            entry["world_size"] for entry in ranks.values()
            if "world_size" in entry
        }
        if len(world_sizes) != 1:
            raise ValueError(f"inconsistent or missing MPI world size in {log}")
        world_size = world_sizes.pop()
        if len(complete) != world_size:
            raise ValueError(
                f"{log}: found hashes for {len(complete)}/{world_size} ranks"
            )
        required_samples = (
            "stage1_seconds",
            "stage2_seconds",
            "cpe_part3_seconds",
            "fread_seconds",
            "fread_mib_per_second",
            "claimed_chunks",
            "scheduler_rma_seconds",
            "scheduler_stage2_seconds",
            "slowest_chunks",
        )
        for key in required_samples:
            if len(samples[key]) != world_size:
                raise ValueError(
                    f"{log}: parsed {len(samples[key])}/{world_size} "
                    f"samples for {key}"
                )
        aggregate_sum = sum(entry["hash_sum"] for entry in complete) & ((1 << 64) - 1)
        aggregate_xor = 0
        for entry in complete:
            aggregate_xor ^= entry["hash_xor"]
        stage1_values = samples["stage1_seconds"]
        stage2_values = samples["stage2_seconds"]
        part3_values = samples["cpe_part3_seconds"]
        rma_values = samples["scheduler_rma_seconds"]
        bandwidth_values = samples["fread_mib_per_second"]
        claimed_values = samples["claimed_chunks"]
        claimed_chunk_total = sum(claimed_values)
        chunk_details_complete = len(chunks) == claimed_chunk_total
        slowest_chunk = max(
            samples["slowest_chunks"], key=lambda item: item[1], default=None
        )
        mapped_stage2 = sum(
            "stage2_seconds" in entry for entry in ranks.values()
        )
        mapped_rma = sum(
            "scheduler_rma_seconds" in entry for entry in ranks.values()
        )
        if stage2_values and min(stage2_values) > 0.0:
            rma_upper_bounds.append(max(rma_values) / min(stage2_values))
        run_rows.append(
            {
                "log": log.name,
                "dataset": dataset,
                "read_mode": read_mode,
                "repeat": repeat,
                "ranks": world_size,
                "calls": sum(entry.get("hash_calls", 0) for entry in complete),
                "bytes": sum(entry.get("hash_bytes", 0) for entry in complete),
                "hash_sum": f"{aggregate_sum:016x}",
                "hash_xor": f"{aggregate_xor:016x}",
                "hash_enabled": (
                    int(all(entry.get("hash_enabled", 1) for entry in complete))
                    if complete else ""
                ),
                "claimed_chunks_range": (
                    f"{min(claimed_values)}-{max(claimed_values)}"
                    if claimed_values else "-"
                ),
                "stage1_range": range_text(stage1_values),
                "stage2_range": range_text(stage2_values),
                "stage2_max_min": f"{ratio(stage2_values):.3f}",
                "part3_range": range_text(part3_values),
                "part3_max_min": f"{ratio(part3_values):.3f}",
                "scheduler_rma_range": range_text(rma_values, 6),
                "fread_bandwidth_range": range_text(bandwidth_values),
                "slowest_chunk": (
                    slowest_chunk[0] if slowest_chunk else ""
                ),
                "slowest_chunk_seconds": (
                    f"{slowest_chunk[1]:.6f}"
                    if slowest_chunk else ""
                ),
                "parsed_chunks": len(chunks),
                "claimed_chunks": claimed_chunk_total,
                "chunk_details_complete": int(chunk_details_complete),
                "rank_details_complete": int(
                    chunk_details_complete and
                    mapped_stage2 == world_size and
                    mapped_rma == world_size
                ),
            }
        )
        chunk_rows.extend(chunks)

    rank_fields = list(rank_rows[0]) if rank_rows else []
    chunk_fields = list(chunk_rows[0]) if chunk_rows else []
    run_fields = list(run_rows[0]) if run_rows else []
    if rank_rows:
        write_tsv(args.result_root / "rank_summary.tsv", rank_rows, rank_fields)
    if chunk_rows:
        write_tsv(args.result_root / "chunk_summary.tsv", chunk_rows, chunk_fields)
    if run_rows:
        write_tsv(args.result_root / "run_summary.tsv", run_rows, run_fields)

    report = args.result_root / "README.md"
    with report.open("w", encoding="utf-8") as handle:
        handle.write("# MPI discard stage2 profile\n\n")
        handle.write(
            "Fixed configuration: `MPI_EXACT_READ_INDEX=1`, "
            "`MPI_INPUT_MODE=dynamic`, `OUTPUT_MODE=discard`, `-1`, `-v 4`. "
            "Hashes are combined offline across ranks with modulo-2^64 sum "
            "and XOR.\n\n"
        )
        handle.write("## Measurement semantics\n\n")
        handle.write(
            "- `stage2` and `CPE part3` are per-rank accumulated active time; "
            "they do not include time for another rank to finish.\n"
            "- Stage1 bandwidth covers timed pipeline `fread` calls. The exact "
            "read-index scan and reference loading happen before the pipeline.\n"
            "- Discard mode creates no SAM file and no output RMA window. With "
            "hashing enabled, each non-empty SAM blob passed to "
            "`swbwa_output_write` is hashed once.\n"
            "- Per-rank hashes are printed independently. The checker combines "
            "them without MPI, using addition modulo 2^64 and XOR; both are "
            "independent of rank assignment and write order.\n"
            "- With `-1`, stage3 hashing is outside stage2, although its cost "
            "still occurs before that rank asks the scheduler for its next "
            "chunk.\n"
            "- Run-level ranges are collected independently of rank labels, "
            "so they remain valid when a batch system interleaves stderr from "
            "different nodes. `rank_summary.tsv` and `chunk_summary.tsv` are "
            "marked incomplete when that interleaving truncates or reorders "
            "detail rows.\n\n"
        )
        handle.write("## Run summary\n\n")
        handle.write(
            "| Dataset | Mode | Run | Chunks/rank | Stage1 (s) | "
            "Stage2 (s) | max/min | CPE part3 (s) | RMA (s) | "
            "Slowest chunk | Chunk rows | fread MiB/s | Hash sum | Hash XOR |\n"
        )
        handle.write(
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|\n"
        )
        for row in run_rows:
            handle.write(
                f"| {row['dataset']} | {row['read_mode']} | {row['repeat']} | "
                f"{row['claimed_chunks_range']} | {row['stage1_range']} | "
                f"{row['stage2_range']} | {row['stage2_max_min']} | "
                f"{row['part3_range']} | {row['scheduler_rma_range']} | "
                f"{row['slowest_chunk']} / {row['slowest_chunk_seconds']} s | "
                f"{row['parsed_chunks']}/{row['claimed_chunks']} | "
                f"{row['fread_bandwidth_range']} | `{row['hash_sum']}` | "
                f"`{row['hash_xor']}` |\n"
            )

        handle.write("\n## Main observations\n\n")
        groups = {}
        for row in run_rows:
            groups.setdefault((row["dataset"], row["read_mode"]), []).append(row)
        matching_groups = 0
        repeated_groups = 0
        for rows in groups.values():
            if len(rows) < 2 or not all(row["hash_enabled"] for row in rows):
                continue
            repeated_groups += 1
            fingerprints = {
                (row["calls"], row["bytes"], row["hash_sum"], row["hash_xor"])
                for row in rows
            }
            if len(fingerprints) == 1:
                matching_groups += 1
        if repeated_groups:
            handle.write(
                f"- Repeated hash-enabled cases with matching aggregate "
                f"fingerprints: {matching_groups}/{repeated_groups}.\n"
            )

        if rma_upper_bounds:
            handle.write(
                f"- Largest scheduler-RMA/stage2 upper bound: "
                f"{max(rma_upper_bounds) * 100.0:.2f}% "
                f"(max RMA divided by min stage2 within a run).\n"
            )

        runs_with_incomplete_details = [
            row for row in run_rows
            if not row["chunk_details_complete"] or
               not row["rank_details_complete"]
        ]
        if runs_with_incomplete_details:
            handle.write(
                f"- Runs with incomplete rank/chunk detail parsing due to "
                f"combined-log interleaving: "
                f"{len(runs_with_incomplete_details)}/{len(run_rows)}. "
                f"Run-level ranges and hashes remain complete.\n"
            )

        if run_rows:
            slowest = max(
                run_rows, key=lambda item: float(item["slowest_chunk_seconds"])
            )
            handle.write(
                f"- Slowest observed chunk: {slowest['dataset']} "
                f"{slowest['read_mode']} chunk {slowest['slowest_chunk']} at "
                f"{slowest['slowest_chunk_seconds']} s.\n"
            )
        bandwidths = [
            row["fread_mib_per_second"] for row in rank_rows
            if row["fread_mib_per_second"] != ""
        ]
        if bandwidths:
            handle.write(
                f"- Lowest measured pipeline fread bandwidth: "
                f"{min(bandwidths):.3f} MiB/s.\n"
            )
        handle.write(
            "\nDetailed data are in `rank_summary.tsv` and "
            "`chunk_summary.tsv`; the original per-rank reports remain in "
            "the `.runNN.log` files.\n"
        )

    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
