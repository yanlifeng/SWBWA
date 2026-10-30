#!/usr/bin/env python3
"""Real gzip reader + processor handoff + pipeline queues + ordered MPI writer.

The processor generates deterministic unmapped SAM, not real alignments.
Both plain and real gzip output exercise the shared processor handoff.
Prepared-output stress cases free SAM before the writer consumes the batch.
All jobs run locally with the fork launcher; fixtures and SAMs are temporary.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile

from test_mpi_gzip_input import build as build_reader
from test_sam_stream_md5 import build_md5_object

ROOT = Path(__file__).resolve().parents[1]


def build(work, mode, sanitizer=None):
    build_reader(work, sanitizer, mode)
    san = (["-fsanitize=" + sanitizer, "-fno-sanitize-recover=all",
            "-fno-omit-frame-pointer", "-fno-pie"] if sanitizer else [])
    flags = ["-O2", "-g", "-pthread", "-I" + str(ROOT / "include"),
             "-I" + str(ROOT / "src/gzip"), "-DSWBWA_USE_MPI=1",
             "-DSWBWA_ENABLE_GZIP_INPUT=1", "-DSWBWA_ENABLE_GZIP_OUTPUT=1",
             "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
             "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_" + mode.upper(),
             "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_ORDERED",
             "-include", str(ROOT / "include/swbwa_config.h"), *san]
    objects = [str(work / (s + ".o")) for s in (
        "mpi_input", "mpi_decoder", "fastq_frame", "codec_job", "backend_host",
        "rapid_firstpass", "swbwa_host_workers")]
    for source in ("src/host/kthread.cpp", "src/host/swbwa_mpi.c", "src/host/swbwa_output.c",
                   "src/gzip/output.c",
                   "tests/test_mpi_gzip_dynamic_pipeline.cpp"):
        src = ROOT / source
        obj = work / (src.stem + ".o")
        cpp = src.suffix == ".cpp"
        compiler = os.getenv("MPICXX", "mpicxx") if cpp else os.getenv("MPICC", "mpicc")
        subprocess.run([compiler, "-std=c++11" if cpp else "-std=gnu11", *flags,
                        "-c", str(src), "-o", str(obj)], check=True)
        objects.append(str(obj))
    binary = work / "test-dynamic-pipeline"
    subprocess.run([os.getenv("MPICXX", "mpicxx"), *objects, *build_md5_object(work),
                    *san, *(["-no-pie"] if sanitizer else []),
                    "-Wl,--wrap=swbwa_mpi_decoder_config_init", "-Wl,--wrap=swbwa_codec_launch",
                    "-Wl,--wrap=swbwa_mpi_decoder_next",
                    "-Wl,--wrap=MPI_Fetch_and_op", "-Wl,--wrap=MPI_Win_flush", "-Wl,--wrap=MPI_Alltoallv",
                    "-Wl,--wrap=pwritev",
                    "-Wl,--wrap=pthread_key_create", "-Wl,--wrap=pthread_setspecific",
                    "-pthread", "-o", str(binary)], check=True)
    return binary


def reports(stderr, prefix):
    rows = re.findall(r"\[" + re.escape(prefix) + r"\] ([^\n]+)", stderr)
    return [dict(field.split("=", 1) for field in row.split() if "=" in field) for row in rows]


def check_stress(case):
    scheduler = re.findall(r"\[MPI gzip rank=(\d+) scheduler\] ([^\n]+)", case["stderr"])
    metrics = {int(rank): dict(field.split("=", 1) for field in row.split())
               for rank, row in scheduler}
    assert len(metrics) == case["ranks"], case
    for rank, row in metrics.items():
        observed = reports(case["stderr"], f"test pipeline rank={rank}")
        assert len(observed) == 1, case
        row = {key: int(value) for key, value in row.items() if key in (
            "direct_claims", "local_claims", "remote_claims", "windows", "data_window_builds",
            "window_capacity_bytes", "fetched_bytes", "decode_drains")}
        count = int(observed[0]["chunks"])
        assert row["local_claims"] + row["remote_claims"] == count, case
        if count:
            assert row["direct_claims"] > 0, case
        if case["stress"] == "single-direct":
            assert row["direct_claims"] == count and row["remote_claims"] == 0, case
            assert all(row[key] == 0 for key in (
                "windows", "data_window_builds", "window_capacity_bytes", "fetched_bytes", "decode_drains")), case
            assert all(float(metrics[rank][key]) == 0 for key in (
                "pack_s", "local_copy_s", "claim_rma_s", "drain_wait_s")), case
            assert int(observed[0]["codec_with_pending"]) > 0, case
        elif case["stress"] in ("decoded-refills", "window-reuse", "window-growth"):
            assert row["windows"] >= 4, case
            if case["stress"] == "window-reuse":
                assert row["data_window_builds"] == 1, case
                assert row["window_capacity_bytes"] >= 1 << 20, case
                assert row["windows"] > row["data_window_builds"], case
    if case["stress"] == "decoded-refills":
        observed = reports(case["stderr"], "test pipeline rank=0")[0]
        assert int(observed["max_overlap_streak"]) >= 2, case
        assert int(observed["decodes_after_overlap"]) > 0, case
        if case["drain"] == 1:
            assert int(metrics[0]["decode_drains"]) > 0, case
    if case["stress"] == "prefix-full":
        first = reports(case["stderr"], "test pipeline rank=0")[0]
        second = reports(case["stderr"], "test pipeline rank=1")[0]
        assert int(first["prefix_releases"]) == 1, case
        assert int(second["prefix_waits"]) > 0 and int(second["decode_full_queue"]) > 0, case
        assert int(second["decode_with_pending"]) > 0 and int(second["codec_with_pending"]) > 0, case
    if case["stress"] in ("exchange-bulk", "exchange-bulk-chunks"):
        observed = [reports(case["stderr"], f"test pipeline rank={rank}")[0]
                    for rank in range(case["ranks"])]
        assert int(observed[0]["bulk_packets"]) >= case["ranks"] - 1, case
        assert int(observed[0]["empty_receives"]) > 0, case
        assert all(int(row["empty_sends"]) > 0 for row in observed[1:]), case
        if case["stress"] == "exchange-bulk-chunks":
            assert int(observed[0]["bulk_payloads"]) >= case["ranks"] - 1, case
    if case["stress"] == "window-reuse":
        assert sum(int(row["remote_claims"]) for row in metrics.values()) > 0, case
    if case["stress"] == "window-growth":
        history = {}
        for rank, fields in scheduler:
            row = dict(field.split("=", 1) for field in fields.split())
            history.setdefault(int(rank), []).append(tuple(int(row[key]) for key in (
                "windows", "data_window_builds", "window_capacity_bytes")))
        for rank, rows in history.items():
            capacity = (2 if rank == 0 else 1) << 20
            expected = ([(0, 0, 0), (1, 1, 1 << 20), (2, 1, 1 << 20), (3, 2, capacity)]
                        if case['early_first'] else
                        [(1, 1, 1 << 20), (2, 1, 1 << 20), (3, 2, capacity), (4, 2, capacity)])
            assert rows[:4] == expected, case
            assert int(metrics[rank]["data_window_builds"]) == 2, case
        case["window_history"] = history
    case["scheduler_metrics"] = metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--input-mode", choices=("static", "dynamic"), default="dynamic")
    parser.add_argument("--focused", action="store_true", help="run only the scheduler regressions")
    parser.add_argument("--early-first", action="store_true", help="deliver the reserved chunk before preparing the shared suffix")
    parser.add_argument("--prefix-ahead", action="store_true", help="publish gzip output lengths independently of disk writes")
    parser.add_argument("--prepare-ahead", action="store_true", help="prepare a bounded FASTQ window on a nested producer")
    parser.add_argument("--parallel-history", action="store_true", help="collectively exchange independent decoder tails")
    parser.add_argument("--binary", type=Path, help="reuse an already built integration executable")
    parser.add_argument("--sanitize", choices=("undefined", "address,undefined"))
    args = parser.parse_args()
    if args.focused and args.input_mode != "dynamic":
        parser.error("--focused requires --input-mode dynamic")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve() if args.binary else build(work, args.input_mode, args.sanitize)
    rng = random.Random(509)
    records = []
    for i in range(384):
        pair = []
        for mate in range(2):
            n = 35 + (i * 11 + mate * 7) % 151
            seq = bytes(rng.choices(b"ACGT", k=n))
            qual = bytes(rng.choices(b"!+5:@AIJ", k=n))
            name = f"read{i:08d}".encode()
            pair.append((name, seq, qual))
        records.append(pair)
    logs = []
    with tempfile.TemporaryDirectory(prefix="fixtures-", dir=work) as tmp:
        tmp = Path(tmp)
        cases = [(ranks, ready, pipeline, mates, n, "", 0)
                 for ranks in (1, 2, 3, 6) for pipeline in ("single", "queue")
                 for mates in (1, 2) for ready in (1, 4) for n in (227,)] if not args.focused else []
        if not args.focused:
            cases += [(3, 4, pipeline, 2, n, "", 0) for pipeline in ("single", "queue") for n in (0, 1)]
        cases = [(*case, compressed, drain, None, 0) for case in cases
                 for compressed in (False, True) for drain in (0, 1)]
        if args.input_mode == "dynamic":
            # Expand only scheduler regressions, not the full ordinary input matrix.
            scheduler_cases = [(2, 1, "queue", 2, 192, "decoded-refills", 64),
                               (3, 2, "queue", 2, 384, "window-reuse", 128),
                               (1, 4, "queue", 2, 24, "single-direct", 8),
                               (2, 2, "queue", 2, 16, "window-growth", 0)]
            cases += [(*case, compressed, drain, None, prepare) for case in scheduler_cases
                      for compressed, prepare in ((False, 0), (True, 0), (True, 1))
                      for drain in (0, 1)]
            cases += [(2, 4, "single", 2, 1, "pending-close", 0, compressed, drain, None, 0)
                      for compressed in (False, True) for drain in (0, 1)]
            cases += [(2, 1, "queue", 2, 32, "prefix-full", 2, True, drain, None, prepare)
                      for drain in (None, 0) for prepare in (0, 1)]
            cases += [(2, 1, "queue", 2, n, "eof-codec", 0, True, drain, None, prepare)
                      for n in (0, 1) for drain in (0, 1) for prepare in (0, 1)]
            cases += [(ranks, 1, "queue", 2, 8, stress, 0, True, 0, None, prepare)
                      for stress in ("exchange-bulk", "exchange-bulk-chunks")
                      for ranks in (2, 3) for prepare in (0, 1)]
        cases += [(2, 1, "single", 2, 8, "", 0, True, drain, None, prepare)
                  for drain in (0, 1) for prepare in (0, 1)]
        cases += [(2, 1, pipeline, 2, 8, "", 0, False, 0, None, 1)
                  for pipeline in ("single", "queue")]
        cases += [(ranks, 1, pipeline, 2, 32, "", 0, compressed, 0, None, None)
                  for ranks in (1, 3) for pipeline in ("single", "queue")
                  for compressed in (False, True)]
        cases += [(2, 4, "queue", 2, 1, "", 0, True, value,
                   "SWBWA_GZIP_DECODE_DRAIN must be 0 or 1", 0)
                  for value in ("", "2", "-1", "true", "00", " 0")]
        cases += [(2, 4, "queue", 2, 1, "drain-mismatch", 0, True, 0,
                   "MPI gzip reader configuration mismatch", 0)]
        cases += [(2, 4, "queue", 2, 1, "early-first-mismatch", 0, True, 0,
                   "MPI gzip reader configuration mismatch", 0)]
        cases += [(2, 4, "queue", 2, 1, "prepare-ahead-mismatch", 0, True, 0,
                   "MPI gzip reader configuration mismatch", 0)]
        if args.prefix_ahead:
            # prefix-full injects a Writer-thread RMA stall. The new publisher
            # deliberately removes that dependency; test blocked pwrite instead.
            cases = [case for case in cases if case[5] != 'prefix-full']
            if args.input_mode == 'dynamic':
                cases += [(2, 1, 'queue', 2, 8, 'prefix-ahead', 0, True, 0, None, 1)]
        for ranks, ready, pipeline, mates, n, stress, member_records, compressed, drain, error, prepare in cases:
            selected = records[:n]
            if stress == "window-growth":
                # With one record per chunk, only rank 0's third-window suffix grows.
                selected[10] = [(name, b"A" * 300000, b"I" * 300000)
                                for name, _, _ in selected[10]]
            if stress in ("exchange-bulk", "exchange-bulk-chunks"):
                selected = [[(name, b"A" * 300000, b"I" * 300000) for name, _, _ in pair]
                            for pair in selected]
            paths = [tmp / "R1.gz", tmp / "R2.gz"]
            for m, path in enumerate(paths):
                raw = [b"@" + pair[m][0] + b" comment\n" + pair[m][1] +
                       b"\n+\n" + pair[m][2] + b"\n" for pair in selected]
                # Member boundaries force later decoder epochs without large fixtures.
                stride = member_records or max(1, n)
                path.write_bytes(b"".join(gzip.compress(b"".join(raw[start:start + stride]), mtime=0)
                                          for start in range(0, max(1, n), stride)))
            expected = b"".join(pair[m][0] + b"\t4\t*\t0\t0\t*\t*\t0\t0\t" + pair[m][1] +
                                b"\t" + pair[m][2] + f"\tGI:i:{i * mates + m}\n".encode()
                                for i, pair in enumerate(selected) for m in range(mates))
            output = tmp / ("output.sam.gz" if compressed else "output.sam")
            env = dict(os.environ, SWBWA_GZIP_READY_CHUNKS=str(ready), HYDRA_LAUNCHER="fork",
                       SWBWA_GZIP_PARALLEL_WRITE="1", SWBWA_GZIP_EARLY_FIRST=str(int(args.early_first)))
            env['SWBWA_GZIP_PREFIX_AHEAD'] = str(int(args.prefix_ahead))
            env['SWBWA_GZIP_PREPARE_AHEAD'] = str(int(args.prepare_ahead))
            env['SWBWA_GZIP_PARALLEL_HISTORY'] = str(int(args.parallel_history))
            if stress == 'prefix-ahead': env['SWBWA_GZIP_PARALLEL_WRITE'] = '0'
            if prepare is None:
                env.pop("SWBWA_GZIP_PREPARE_OUTPUT", None)
            else:
                env["SWBWA_GZIP_PREPARE_OUTPUT"] = str(prepare)
            env.setdefault("ASAN_OPTIONS", "detect_leaks=0:halt_on_error=1")
            env["MALLOC_PERTURB_"] = "165"
            env.pop("SWBWA_OUTPUT_MD5", None)
            env.pop("HYDRA_HOST_FILE", None)
            if drain is None:
                env.pop("SWBWA_GZIP_DECODE_DRAIN", None)
            else:
                env["SWBWA_GZIP_DECODE_DRAIN"] = str(drain)
            env["PIPELINE_TEST_STRESS"] = stress
            command = ["mpirun", "-launcher", "fork", "-n", str(ranks), str(binary),
                       str(paths[0]), str(paths[1]) if mates == 2 else "-", str(output), pipeline]
            case = dict(ranks=ranks, ready=ready, pipeline=pipeline, mates=mates, records=n,
                        stress=stress, member_records=member_records, gzip_output=compressed,
                        drain=drain, prepare=prepare, early_first=args.early_first, prefix_ahead=args.prefix_ahead,
                        prepare_ahead=args.prepare_ahead, sanitizer=args.sanitize, command=command)
            try:
                process = subprocess.run(command, env=env, text=True, capture_output=True, timeout=90)
            except subprocess.TimeoutExpired as failure:
                case.update(status="TIMEOUT", stdout=(failure.stdout or b"").decode(errors="replace"),
                            stderr=(failure.stderr or b"").decode(errors="replace"))
                logs.append(case)
                (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
                raise
            case.update(returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
            logs.append(case)
            (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
            if error or stress == "pending-close":
                assert process.returncode != 0, case
                if error:
                    assert error in process.stderr, case
                else:
                    assert "MPI gzip close requires completed output batches" in process.stderr, case
                    assert reports(process.stderr, "test pending close rank=0") == [
                        {"validated_EOF": "1", "pending": "1"}], case
                    assert "pending close unexpectedly returned" not in process.stderr, case
                case["status"] = "EXPECTED_REJECTION"
                (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
                print("EXPECTED_REJECTION", stress or "invalid-drain", repr(drain), compressed, flush=True)
                continue
            assert process.returncode == 0, case
            actual = output.read_bytes()
            assert (gzip.decompress(actual) if compressed else actual) == expected, case
            for rank in range(ranks):
                keys = reports(process.stderr, f"test pipeline keys rank={rank}")
                expected_pipelines = 2 + (pipeline == "queue")
                assert len(keys) == 1 and int(keys[0]["pipelines"]) == expected_pipelines, case
                expected_producers = int(args.prepare_ahead and ranks > 1 and
                                         args.input_mode == 'dynamic' and pipeline == 'queue')
                preparation = reports(process.stderr, f"MPI gzip rank={rank} prepare_ahead")[-1]
                assert int(preparation['enabled']) == int(args.prepare_ahead), case
                assert int(preparation['used']) == expected_producers, case
                if not expected_producers:
                    assert int(preparation['prepared_chunks']) == 0, case
                assert all(int(keys[0][key]) == 2 * expected_pipelines + expected_producers for key in (
                    "sets", "clears", "clear_direct_calls")), case
                observed = reports(process.stderr, f"test pipeline rank={rank}")[0]
                scheduler = reports(process.stderr, f"MPI gzip rank={rank} scheduler")[-1]
                local = reports(process.stderr, f"MPI gzip rank={rank} reader_local")[-1]
                assert local["exchange_pack_mode"] == "plan_and_local_copy", case
                assert local["exchange_prepare_mode"] == "direct_wire_fill", case
                exchange_parts = [float(local[key]) for key in (
                    "exchange_pack_s", "exchange_prepare_s", "exchange_receive_s")]
                assert all(value >= 0 for value in exchange_parts), case
                assert abs(sum(exchange_parts) - float(local["exchange_local_s"])) <= 0.000003, case
                if ranks == 1:
                    assert exchange_parts[1:] == [0, 0], case
                assert int(scheduler["decode_drain_enabled"]) == int(drain or 0), case
                assert int(scheduler['early_first_enabled']) == int(args.early_first), case
                if args.early_first and ranks > 1 and args.input_mode == 'dynamic' and int(observed["chunks"]):
                    assert int(scheduler['early_first_chunks']) > 0, case
                if not drain:
                    assert int(scheduler["decode_drains"]) == 0 and float(scheduler["drain_wait_s"]) == 0, case
                elif args.input_mode == "dynamic" and ranks > 1:
                    assert int(observed["decode_with_pending"]) == int(observed["codec_with_pending"]) == 0, case
                eof = int(stress == "eof-codec")
                assert int(observed["eof_codecs"]) == eof, case
                assert int(observed["output_codecs"]) == (int(observed["chunks"]) if compressed else 0) + eof, case
                prepared = int(observed["chunks"]) if compressed and prepare != 0 else 0
                assert all(int(observed[key]) == prepared for key in (
                    "prepared_chunks", "prepared_writes", "prepare_codecs", "freed_sam_chunks")), case
                if args.prefix_ahead and compressed and prepare != 0 and args.input_mode == 'dynamic':
                    prefix = reports(process.stderr, f'SWBWA gzip prefix rank {rank:06d}/{ranks:06d}')
                    assert len(prefix) == 1 and int(prefix[0]['submitted']) == prepared and int(prefix[0]['error']) == 0, case
            if stress:
                check_stress(case)
            case["status"] = "PASS"
            (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
            output.unlink()
            print("PASS", ranks, ready, pipeline, mates, n, stress, "gzip", compressed,
                  "drain", drain, "prepare", prepare, flush=True)
    print("MPI_GZIP_ORDERED_PIPELINE_PASS", len(logs), flush=True)


if __name__ == "__main__":
    main()
