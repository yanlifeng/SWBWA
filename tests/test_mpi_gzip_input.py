#!/usr/bin/env python3
"""Local real-decoder/frame/C-reader tests. mpirun uses only the fork launcher."""
import argparse
import json
import os
from pathlib import Path
import random
import re
import shlex
import subprocess
import tempfile
import zlib

from test_mpi_gzip_decoder import encoded, stored

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "src" / "gzip"


def add_transfer_argument(parser):
    parser.add_argument("--transfer", "--dynamic-transfer", "--transfer-mode", dest="dynamic_transfer",
                        choices=("staged", "direct", "local"),
                        default=os.environ.get("SWBWA_GZIP_DYNAMIC_TRANSFER", "staged"),
                        help="dynamic payload strategy (default: SWBWA_GZIP_DYNAMIC_TRANSFER or staged)")


def check_chunk_trace(stderr, records, ranks, input_mode):
    estimates = set(map(int, re.findall(r"records_per_chunk=(\d+)", stderr)))
    assert len(estimates) == 1, stderr
    k = estimates.pop()
    assert k > 0, stderr
    previous = [-1] * ranks
    chunks = {}
    for fields in re.findall(
        r"\[MPI gzip rank=(\d+)\] chunk=(\d+) first_record=(\d+) "
        r"end_record_exclusive=(\d+) records=(\d+) R1_bytes=(\d+) R2_bytes=(\d+)", stderr
    ):
        rank, chunk, first, end, count, bytes1, bytes2 = map(int, fields)
        assert 0 <= rank < ranks and chunk > previous[rank], fields
        previous[rank] = chunk
        assert first == min(chunk * k, records) and end == first + count, fields
        assert count == min(k, records - first), fields
        if input_mode == "static":
            assert chunk % ranks == rank, fields
        else:
            assert count > 0, fields
        if count:
            assert chunk not in chunks, fields
            chunks[chunk] = (chunk, first, count, bytes1, bytes2)
    assert sorted(chunks) == list(range((records + k - 1) // k)), chunks
    return [chunks[i] for i in sorted(chunks)]


def check_payload_traffic(stderr, ranks, payload_bytes, direct_source):
    rows = re.findall(r"\[MPI gzip rank=(\d+) reader_local\] ([^\n]+)", stderr)
    stats = {int(rank): dict(re.findall(r"(\w+)=(\S+)", row)) for rank, row in rows}
    assert set(stats) == set(range(ranks)), stderr
    local = sum(int(row["local_payload_bytes"]) for row in stats.values())
    remote = sum(int(row["remote_payload_bytes"]) for row in stats.values())
    assert local + remote == payload_bytes, (local, remote, payload_bytes)
    if direct_source:
        # This counter measures pre-distribution, not on-demand source RMA reads.
        assert all(int(row["remote_payload_bytes"]) == 0 for row in stats.values()), stats
    return remote


def check_source_traffic(stderr, ranks, payload_bytes, transfer):
    rows = re.findall(r"\[MPI gzip rank=(\d+) direct_source\] ([^\n]+)", stderr)
    stats = {int(rank): dict(re.findall(r"(\w+)=(\S+)", row)) for rank, row in rows}
    assert set(stats) == set(range(ranks)), stderr
    enabled = int(ranks > 1 and transfer != "staged")
    for row in stats.values():
        assert row["transfer"] == transfer and int(row["enabled"]) == enabled, row
        assert all(int(row[key]) >= 0 for key in (
            "local_fastq_bytes", "remote_fastq_bytes", "metadata_contribution_bytes", "get_calls",
            "frontier_chunks", "short_windows")), row
        if "flush_calls" in row:
            assert 0 <= int(row["flush_calls"]) <= int(row["get_calls"]), row
            assert bool(int(row["flush_calls"])) == bool(int(row["get_calls"])), row
    if enabled:
        assert sum(int(row["local_fastq_bytes"]) + int(row["remote_fastq_bytes"])
                   for row in stats.values()) == payload_bytes, stats
        for row in stats.values():
            assert bool(int(row["remote_fastq_bytes"])) == bool(int(row["get_calls"])), row
    else:
        assert all(int(row[key]) == 0 for row in stats.values() for key in (
            "local_fastq_bytes", "remote_fastq_bytes", "get_calls", "frontier_chunks", "short_windows")), stats
    return stats


def check_source_frontier(stderr, chunk_count, ranks, ready, transfer, stats, early_first=False):
    if ranks == 1 or transfer == "staged":
        return []
    owners = {int(chunk): int(rank) for rank, chunk in re.findall(
        r"\[MPI gzip rank=(\d+)\] chunk=(\d+)", stderr)}
    windows = [[] for _ in range(ranks)]
    for fields in re.findall(
        r"\[MPI gzip rank=(\d+) source_window\] begin=(\d+) end=(\d+) first=(-?\d+)"
        r"(?: frontier=(\d+))?", stderr
    ):
        rank, begin, end, first = map(int, fields[:4])
        frontier = int(fields[4]) if fields[4] else min(ranks, end - begin)
        assert 0 <= rank < ranks, fields
        windows[rank].append((begin, end, first, frontier))
    boundaries = [(begin, end) for begin, end, _, _ in windows[0]]
    descriptors = [(begin, end, frontier) for begin, end, _, frontier in windows[0]]
    assert all([(begin, end, frontier) for begin, end, _, frontier in row] == descriptors
               for row in windows), windows
    frontier_counts = [0] * ranks
    cursor = cycle_end = short_candidates = 0
    for index, (begin, end, frontier) in enumerate(descriptors):
        reserve = not early_first or begin >= cycle_end
        if early_first and reserve:
            cycle_end = begin + ranks * ready
        limit = (min(cycle_end, begin + ranks) if reserve else cycle_end) if early_first else begin + ranks * ready
        assert begin == cursor and begin < end <= min(limit, chunk_count), descriptors
        assert frontier == (min(ranks, end - begin) if reserve else 0), descriptors
        short_candidates += end < limit
        cursor = end
        reserved = [windows[rank][index][2] for rank in range(ranks)]
        assert all(chunk >= -1 for chunk in reserved), reserved
        assert sorted(chunk for chunk in reserved if chunk >= 0) == list(range(begin, begin + frontier)), reserved
        for rank, chunk in enumerate(reserved):
            if chunk == -1:
                continue
            assert begin <= chunk < end and owners[chunk] == rank, (begin, end, rank, chunk)
            if transfer == "direct":
                assert rank == chunk - begin, (chunk, rank)
            frontier_counts[rank] += 1
    assert cursor == chunk_count, (cursor, chunk_count)
    assert [int(stats[rank]["frontier_chunks"]) for rank in range(ranks)] == frontier_counts, stats
    short = {int(row["short_windows"]) for row in stats.values()}
    assert len(short) == 1 and short.pop() <= short_candidates, stats
    return boundaries


def build(work, sanitizer, input_mode="static"):
    flags = ["-O1" if sanitizer else "-O2", "-g", "-Wall", "-Wextra",
             "-D_FILE_OFFSET_BITS=64", "-D_POSIX_C_SOURCE=200809L",
             "-DSWBWA_CODEC_RAPID=1", "-DSWBWA_CODEC_LDM=2", "-DSWBWA_USE_MPI=1",
             "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_" + input_mode.upper(),
             "-I" + str(ROOT / "include"), "-I" + str(ENGINE),
             "-I" + str(ENGINE / "vendor" / "rapidgzip")]
    san = ["-fsanitize=" + sanitizer, "-fno-sanitize-recover=all", "-fno-omit-frame-pointer", "-fno-pie"] if sanitizer else []
    link_san = [*san, *(["-no-pie"] if sanitizer else [])]
    flags += san
    sources = [ENGINE / name for name in ("mpi_input.cpp", "fastq_frame.cpp", "mpi_decoder.cpp",
               "codec_job.c", "backend_host.c", "rapid_firstpass.cpp")]
    sources += [ROOT / "tests" / "test_mpi_gzip_input.cpp"]
    sources += [ROOT / "src" / "host" / "swbwa_host_workers.c"]
    objects = []
    for source in sources:
        obj = work / (source.stem + ".o")
        cpp = source.suffix == ".cpp"
        compiler = os.environ.get("MPICXX", "mpicxx") if cpp else os.environ.get("CC", "cc")
        standard = "c++17" if source.name == "rapid_firstpass.cpp" else "c++11" if cpp else "c11"
        subprocess.run([compiler, "-std=" + standard, *flags, "-c", str(source), "-o", str(obj)], check=True)
        objects.append(str(obj))
    binary = work / "test-mpi-gzip-input"
    subprocess.run([os.environ.get("MPICXX", "mpicxx"), *objects, *link_san,
                    "-Wl,--wrap=swbwa_mpi_decoder_config_init", "-Wl,--wrap=swbwa_mpi_decoder_next",
                    "-Wl,--wrap=MPI_Alltoallv",
                    "-lz", "-pthread", "-o", str(binary)], check=True)
    frame = work / "test-fastq-frame"
    subprocess.run([os.environ.get("CXX", "c++"), "-std=c++11", *flags,
                    str(ROOT / "tests" / "test_fastq_frame.cpp"), str(work / "fastq_frame.o"),
                    *(["-no-pie"] if sanitizer else []), "-o", str(frame)], check=True)
    subprocess.run([os.environ.get("CC", "cc"), "-std=c11", "-Wall", "-Wextra", "-Werror",
                    "-I" + str(ROOT / "include"), "-x", "c", "-fsyntax-only", "-"],
                   input=b'#include "swbwa_mpi_gzip.h"\n', check=True)
    return binary, frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--sanitize", choices=("undefined", "address,undefined"))
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--input-mode", choices=("static", "dynamic"), default="static")
    add_transfer_argument(parser)
    args = parser.parse_args()
    if args.dynamic_transfer not in ("staged", "direct", "local"):
        parser.error("SWBWA_GZIP_DYNAMIC_TRANSFER must be staged, direct or local")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.binary:
        binary = args.binary.resolve()
        frame = binary.with_name("test-fastq-frame")
    else:
        binary, frame = build(work, args.sanitize, args.input_mode)
    subprocess.run([str(frame)], check=True)
    rng = random.Random(2817)
    reports, logs = [], []
    with tempfile.TemporaryDirectory(prefix="fixtures-", dir=work) as temporary:
        root = Path(temporary)

        def run(name, r1, r2=None, blob1=None, blob2=None, ranks=2,
                target=32768, mode="normal", error=None, defaults=False, lookahead=1 << 20, ready=4,
                workers=384, transfer=None, env_overrides=None):
            transfer = transfer or args.dynamic_transfer
            direct_source = args.input_mode == "dynamic" and transfer != "staged"
            raw1, gz1 = root / (name + ".R1.raw"), root / (name + ".R1.gz")
            raw1.write_bytes(r1); gz1.write_bytes(encoded(r1) if blob1 is None else blob1)
            raw2 = gz2 = "-"
            if r2 is not None:
                raw2, gz2 = root / (name + ".R2.raw"), root / (name + ".R2.gz")
                raw2.write_bytes(r2); gz2.write_bytes(encoded(r2) if blob2 is None else blob2)
            command = ["mpirun", "-launcher", "fork", "-n", str(ranks), str(binary),
                       str(gz1), str(gz2), str(raw1), str(raw2), str(target), str(workers), "1", mode]
            env = dict(os.environ)
            env["SWBWA_GZIP_DYNAMIC_TRANSFER"] = transfer
            if direct_source:
                env.pop("SWBWA_GZIP_EARLY_FIRST", None)
                env.pop("SWBWA_GZIP_PREPARE_AHEAD", None)
            env["FASTQ_TEST_LOOKAHEAD"] = str(lookahead)
            env["SWBWA_GZIP_READY_CHUNKS"] = str(ready)
            env["MALLOC_PERTURB_"] = "165"
            if defaults:
                env["FASTQ_TEST_DEFAULT_CONFIG"] = "1"
            else:
                env.pop("FASTQ_TEST_DEFAULT_CONFIG", None)
            env.setdefault("ASAN_OPTIONS", "detect_leaks=0")
            for key, value in (env_overrides or {}).items():
                if value is None:
                    env.pop(key, None)
                else:
                    env[key] = value
            process = subprocess.run(command, text=True, capture_output=True, timeout=120, env=env)
            logs.append(dict(case=name, transfer=transfer, env_overrides=env_overrides,
                             command=shlex.join(command), returncode=process.returncode,
                             stdout=process.stdout, stderr=process.stderr))
            if error:
                assert process.returncode != 0 and error in process.stderr, logs[-1]
                assert "MPI_FASTQ_RESULT" not in process.stdout, logs[-1]
                reports.append(dict(case=name, status="EXPECTED_REJECTION"))
                print("EXPECTED_REJECTION", name, flush=True)
                return None
            assert process.returncode == 0, logs[-1]
            result_lines = [line[len("MPI_FASTQ_RESULT "):] for line in process.stdout.splitlines()
                            if line.startswith("MPI_FASTQ_RESULT ")]
            assert len(result_lines) == 1, logs[-1]
            result = json.loads(result_lines[0])
            assert result["records"] == r1.count(b"\n") // 4
            if defaults:
                widths = list(map(int, re.findall(r"initial_batch_jobs=(\d+)", process.stderr)))
                assert len(widths) == ranks * (1 + (r2 is not None)), logs[-1]
                assert all(width == workers for width in widths), logs[-1]
            result["chunk_layout"] = check_chunk_trace(process.stderr, result["records"], ranks, args.input_mode)
            result["remote_payload_bytes"] = check_payload_traffic(
                process.stderr, ranks, len(r1) + len(r2 or b""), direct_source)
            exchanges = re.findall(r"\[test consumer rank=(\d+)\] .* alltoallv_calls=(\d+)", process.stderr)
            assert len(exchanges) == ranks and {int(rank) for rank, _ in exchanges} == set(range(ranks)), logs[-1]
            if direct_source:
                assert all(int(count) == 0 for _, count in exchanges), logs[-1]
            result["boundary_records"] = sum(map(int, re.findall(r"boundary_records=(\d+)", process.stderr)))
            result["max_epochs"] = max(map(int, re.findall(r"epochs=(\d+)", process.stderr)))
            result["consumer_chunks"] = dict((int(rank), int(count)) for rank, count in re.findall(
                r"\[test consumer rank=(\d+)\] chunks=(\d+)", process.stderr))
            result["remote_claims"] = sum(map(int, re.findall(r"remote_claims=(\d+)", process.stderr)))
            if args.input_mode == "dynamic":
                result["source_traffic"] = check_source_traffic(
                    process.stderr, ranks, len(r1) + len(r2 or b""), transfer)
                result["source_windows"] = check_source_frontier(
                    process.stderr, len(result["chunk_layout"]), ranks, ready, transfer, result["source_traffic"],
                    early_first=env.get("SWBWA_GZIP_EARLY_FIRST") == "1")
                schedulers = re.findall(r"\[MPI gzip rank=\d+ scheduler\] (.*)", process.stderr)
                assert len(schedulers) == ranks, logs[-1]
                stats = [dict(re.findall(r"(\w+)=(\S+)", line)) for line in schedulers]
                if ranks == 1:
                    assert all(int(s["windows"]) == 0 and int(s["data_window_builds"]) == 0 and
                               int(s["direct_claims"]) == result["consumer_chunks"][0] for s in stats)
                elif transfer == "staged":
                    assert all(int(s["direct_claims"]) <= int(s["local_claims"]) for s in stats)
                    assert all(int(s["direct_claims"]) <= int(s["windows"]) and
                               int(s["data_window_builds"]) <= int(s["windows"]) for s in stats)
                else:
                    assert all(int(s["direct_claims"]) == 0 for s in stats), logs[-1]
                if ready == 1 and transfer == "staged":
                    assert result["remote_claims"] == 0
                    assert all(int(s["data_window_builds"]) == 0 for s in stats)
                result["direct_claims"] = sum(int(s["direct_claims"]) for s in stats)
                result["data_window_builds"] = max(int(s["data_window_builds"]) for s in stats)
            reports.append(dict(case=name, transfer=transfer, status="PASS", result=result))
            print("PASS", name, transfer, {key: value for key, value in result.items()
                  if key not in ("chunk_layout", "source_traffic", "source_windows")}, flush=True)
            return result

        r1parts, r2parts = [], []
        for i in range(9000):
            n = 35 + i * 17 % 150
            seq1 = bytes(rng.choices(b"ACGTN", k=n))
            qual1 = bytes(rng.choices(b"!+5:@AIJ", k=n))
            r1parts.append(b"@read%08d mate1\n" % i + seq1 + b"\n+note\n" + qual1 + b"\n")
            r2parts.append(b"@read%08d mate2\n" % i + b"A" * n + b"\n+note\n" + b"@" * n + b"\n")
        r1, r2 = b"".join(r1parts), b"".join(r2parts)
        tiny = b"@one comment\nAC\n+\n@+\n"
        try:
            run("empty-pair", b"", b"", ranks=3)
            assert run("one-record-empty-ranks", tiny, tiny, ranks=3)["empty_ready"] == (2 if args.input_mode == "static" else 0)
            for ranks in (1, 2, 3):
                result = run(f"pair-skew-ranks{ranks}", r1, r2, ranks=ranks)
                assert result["max_epochs"] > 1 and result["boundary_records"] > 0
            result = run("pair-skew-swapped-ranks3", r2, r1, ranks=3)
            assert result["max_epochs"] > 1 and result["boundary_records"] > 0
            for workers in (64, 384):
                run(f"default-decoder-config-workers{workers}", r1, r2, defaults=True, workers=workers)
            run("thread-multiple-reader", r1, r2, mode="thread")
            run("single-end", r1, ranks=3)
            for name, strategy in (("fixed", zlib.Z_FIXED), ("huffman", zlib.Z_HUFFMAN_ONLY), ("rle", zlib.Z_RLE)):
                run(name, r1, r2, blob1=encoded(r1, strategy), blob2=encoded(r2, strategy), lookahead=4 << 20)
            run("stored", r1, r2, blob1=stored(r1), blob2=stored(r2))
            concat1 = encoded(b"") + encoded(r1[:43]) + encoded(r1[43:451231]) + encoded(b"") + encoded(r1[451231:])
            concat2 = encoded(r2[:10101]) + encoded(b"") + encoded(r2[10101:]) + encoded(b"")
            run("concat-members-cross-record", r1, r2, blob1=concat1, blob2=concat2, ranks=3)
            result = run("concat-members-cross-record-swapped", r2, r1, blob1=concat2, blob2=concat1, ranks=3)
            assert result["boundary_records"] > 0
            run("exact-name-different-comments", b"@same a\nAC\n+\n++\n", b"@same\tb\nAC\n+\n@@\n")
            run("record-chunks-one", b"".join(r1parts[:101]), b"".join(r2parts[:101]), target=1, ranks=3)
            long_name = bytes(rng.choices(b"ACGT", k=100000))
            long_record = b"@" + long_name + b"\nAC\n+\n@@\n"
            result = run("no-record-first-epoch", long_record, long_record,
                         blob1=stored(long_record, 997), blob2=stored(long_record, 3991))
            assert result["max_epochs"] > 1 and result["boundary_records"] == 2
            run("bad-header", b"!one\nAC\n+\nII\n", mode="reject", error="invalid FASTQ")
            run("bad-plus", b"@one\nAC\n-\nII\n", mode="reject", error="invalid FASTQ")
            run("bad-quality-length", b"@one\nAC\n+\nI\n", mode="reject", error="invalid FASTQ")
            run("missing-final-lf", tiny[:-1], mode="reject", error="truncated FASTQ")
            run("crlf", tiny.replace(b"\n", b"\r\n"), mode="reject", error="invalid FASTQ")
            run("nul", tiny.replace(b"AC", b"A\0"), mode="reject", error="invalid FASTQ")
            run("name-mismatch", tiny, tiny.replace(b"one", b"two"), mode="reject", error="name token mismatch")
            run("no-slash-normalization", b"@n/1\nAC\n+\nII\n", b"@n/2\nAC\n+\nII\n",
                mode="reject", error="name token mismatch")
            run("mate-count-mismatch", r1, b"".join(r2parts[:-1]), mode="reject", error="record count mismatch")
            run("empty-mate-mismatch", b"", tiny, mode="reject", error="record count mismatch")
            run("unequal-layout", tiny, b"@one comment\nACG\n+\n@@@\n")
            incompatible = []
            for i, record in enumerate(r2parts):
                lines = record.splitlines()
                delta = 15 if i < 4500 else -15
                n = len(lines[1]) + delta
                incompatible.append(lines[0] + b"\n" + b"A" * n + b"\n" + lines[2] + b"\n" + b"@" * n + b"\n")
            run("equal-bytes-independent-partitions", r1, b"".join(incompatible), target=4 << 20)
            interior = list(r1parts)
            interior[17] = interior[17].replace(b"+note", b"-note")
            run("invalid-interior-record", b"".join(interior), r2, target=4 << 20,
                mode="reject", error="invalid FASTQ")
            interior = list(r2parts)
            interior[17] = interior[17].replace(b"read00000017", b"read00000018")
            run("interior-name-mismatch", r1, b"".join(interior), target=4 << 20,
                mode="reject", error="name token mismatch")
            oversized = b"@x\n" + b"A" * (1 << 19) + b"\n+\n" + b"I" * (1 << 19) + b"\n"
            run("record-cap", oversized, blob1=stored(oversized), mode="reject", error="cap")
            corrupt = bytearray(encoded(r2)); corrupt[-8] ^= 1
            run("mate2-crc", r1, r2, blob2=bytes(corrupt), mode="reject", error="CRC/ISIZE mismatch")
            corrupt = bytearray(encoded(tiny)); corrupt[-8] ^= 1
            run("late-concat-crc", r1 + tiny, blob1=encoded(r1) + corrupt, mode="reject", error="CRC/ISIZE mismatch")
            run("config-mismatch", tiny, tiny, mode="config-mismatch", error="reader configuration mismatch")
            run("early-close", tiny, tiny, mode="early-close", error="requires globally validated EOF")
            if args.input_mode == "dynamic":
                run("default-transfer-staged", tiny, transfer="staged",
                    env_overrides={"SWBWA_GZIP_DYNAMIC_TRANSFER": None})
                baseline = None
                if args.dynamic_transfer != "staged":
                    baseline = run("staged-layout-reference", r1, r2, ranks=3, transfer="staged")
                    assert baseline["remote_payload_bytes"] > 0
                    concat_baseline = run("staged-concat-reference", r2, r1, blob1=concat2,
                                          blob2=concat1, ranks=3, transfer="staged")
                    concat_result = run("source-concat-layout", r2, r1, blob1=concat2,
                                        blob2=concat1, ranks=3)
                    assert concat_result["chunk_layout"] == concat_baseline["chunk_layout"]
                for ready in (1, 2, 4, 8):
                    result = run(f"ready-window-{ready}", r1, r2, ranks=3, ready=ready)
                    if baseline is not None:
                        assert result["chunk_layout"] == baseline["chunk_layout"]
                    if ready == 1 and args.dynamic_transfer == "local":
                        # Every claim is in the bijective frontier, so this
                        # comparison is independent of consumer timing.
                        direct = run("direct-frontier-reference", r1, r2, ranks=3, ready=1, transfer="direct")
                        assert result["chunk_layout"] == direct["chunk_layout"]
                        assert sum(int(row["remote_fastq_bytes"]) for row in result["source_traffic"].values()) < sum(
                            int(row["remote_fastq_bytes"]) for row in direct["source_traffic"].values())
                result = run("skewed-consumer", r1, r2, ranks=3, mode="skew", ready=8)
                if args.dynamic_transfer == "staged":
                    assert result["remote_claims"] > 0
                else:
                    assert result["chunk_layout"] == baseline["chunk_layout"]
                counts = result["consumer_chunks"]
                assert counts[0] + counts[1] > 2 * counts[2], counts
                run("bad-ready-limit", tiny, ready=0, mode="reject", error="SWBWA_GZIP_READY_CHUNKS")
                run("ready-config-mismatch", tiny, mode="window-mismatch", error="reader configuration mismatch")
                run("invalid-completion", tiny, mode="bad-complete", error="completion out of order")
                run("transfer-config-mismatch", tiny, mode="transfer-mismatch", error="reader configuration mismatch")
                for value in ("", "DIRECT", "invalid", "direct "):
                    run(f"invalid-transfer-{value!r}", tiny, mode="reject", error="SWBWA_GZIP_DYNAMIC_TRANSFER",
                        env_overrides={"SWBWA_GZIP_DYNAMIC_TRANSFER": value})
                if args.dynamic_transfer != "staged":
                    early = {"SWBWA_GZIP_EARLY_FIRST": "1"}
                    for ready in (1, 2, 4, 8):
                        result = run(f"early-skew-ready{ready}", r1, r2, blob1=concat1, blob2=concat2,
                                     ranks=3, ready=ready, mode="skew", env_overrides=early)
                        assert result["chunk_layout"] == baseline["chunk_layout"]
                        assert result["consumer_chunks"][2] > 0
                        if ready > 1:
                            assert "frontier=0" in logs[-1]["stderr"], logs[-1]
                        for n in (0, 1, 2):
                            run(f"early-empty-ranks-ready{ready}-records{n}", b"".join(r1parts[:n]),
                                b"".join(r2parts[:n]), ranks=3, ready=ready, target=1, env_overrides=early)
                    # Single final stored blocks avoid a second discovery read
                    # (which needs 64 KiB headroom even for an empty end block).
                    pressure = dict(blob1=b"".join(encoded(r1[i:i + 16381], level=0) for i in range(0, len(r1), 16381)),
                                    blob2=b"".join(encoded(r2[i:i + 24571], level=0) for i in range(0, len(r2), 24571)),
                                    ranks=3, ready=16, target=65536)
                    reference = run("staged-pressure-layout-reference", r1, r2, transfer="staged", **pressure)
                    for early_first in (0, 1):
                        result = run("source-queue-pressure" + ("-early" if early_first else ""), r1, r2,
                                     env_overrides={"FASTQ_TEST_DECODED_BYTES": "65536",
                                                    "SWBWA_GZIP_EARLY_FIRST": str(early_first)}, **pressure)
                        assert result["chunk_layout"] == reference["chunk_layout"]
                        assert all(int(row["short_windows"]) > 0 for row in result["source_traffic"].values())
                        assert any(end - begin < 3 * 16 and end < len(result["chunk_layout"])
                                   for begin, end in result["source_windows"])
                        if early_first:
                            # Pressure must split the ticket-only suffix more than once.
                            assert len(re.findall(r"rank=0 source_window\].* first=-1 frontier=0", logs[-1]["stderr"])) > 1
                        limits = set(map(int, re.findall(r"queue_limit_per_mate_bytes=(\d+)", logs[-1]["stderr"])))
                        peaks = list(map(int, re.findall(r"queue_peak_accounted_bytes=(\d+)", logs[-1]["stderr"])))
                        assert limits == {4 * 65536} and len(peaks) == 6, logs[-1]
                        assert 2 * 65536 < max(peaks) <= 4 * 65536, peaks
                    run("unsupported-SWBWA_GZIP_PREPARE_AHEAD", tiny, mode="reject", error="PREPARE_AHEAD",
                        env_overrides={"SWBWA_GZIP_PREPARE_AHEAD": "1"})
                    for setting in ("SWBWA_GZIP_PREPARE_AHEAD", "SWBWA_GZIP_EARLY_FIRST"):
                        run(f"single-rank-{setting}", tiny, ranks=1, env_overrides={setting: "1"})
        finally:
            (work / "test-results.json").write_text(json.dumps(reports, indent=2) + "\n")
            (work / "test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
    print("MPI_FASTQ_INPUT_TESTS_PASS", len(reports), flush=True)


if __name__ == "__main__":
    main()
