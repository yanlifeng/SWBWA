#!/usr/bin/env python3
"""Real gzip reader + processor handoff + pipeline queues + ordered MPI writer.

The processor generates deterministic unmapped SAM, not real alignments.
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


def build(work, mode):
    build_reader(work, None, mode)
    flags = ["-O2", "-g", "-pthread", "-I" + str(ROOT / "include"),
             "-I" + str(ROOT / "src/gzip"), "-DSWBWA_USE_MPI=1",
             "-DSWBWA_ENABLE_GZIP_INPUT=1", "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
             "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_" + mode.upper(),
             "-DSWBWA_OUTPUT_MODE=SWBWA_OUTPUT_SINGLE_ORDERED",
             "-include", str(ROOT / "include/swbwa_config.h")]
    objects = [str(work / (s + ".o")) for s in (
        "mpi_input", "mpi_decoder", "fastq_frame", "codec_job", "backend_host",
        "rapid_firstpass", "swbwa_host_workers")]
    for source in ("src/host/kthread.cpp", "src/host/swbwa_mpi.c", "src/host/swbwa_output.c",
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
                    "-Wl,--wrap=swbwa_mpi_decoder_config_init", "-Wl,--wrap=swbwa_codec_launch",
                    "-Wl,--wrap=swbwa_mpi_decoder_next",
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
        assert row["direct_claims"] > 0, case
        if case["stress"] == "single-direct":
            assert row["direct_claims"] == count and row["remote_claims"] == 0, case
            assert all(row[key] == 0 for key in (
                "windows", "data_window_builds", "window_capacity_bytes", "fetched_bytes", "decode_drains")), case
            assert all(float(metrics[rank][key]) == 0 for key in (
                "pack_s", "local_copy_s", "claim_rma_s", "drain_wait_s")), case
            assert int(observed[0]["codec_with_pending"]) > 0, case
        else:
            assert row["windows"] >= 4, case
            if case["stress"] == "window-reuse":
                assert row["data_window_builds"] == 1, case
                assert row["window_capacity_bytes"] >= 1 << 20, case
                assert row["windows"] > row["data_window_builds"], case
    if case["stress"] == "decoded-refills":
        observed = reports(case["stderr"], "test pipeline rank=0")[0]
        assert int(observed["max_overlap_streak"]) >= 2, case
        assert int(observed["decodes_after_overlap"]) > 0, case
        assert int(metrics[0]["decode_drains"]) > 0, case
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
            assert rows[:4] == [(1, 1, 1 << 20), (2, 1, 1 << 20),
                                (3, 2, capacity), (4, 2, capacity)], case
            assert int(metrics[rank]["data_window_builds"]) == 2, case
        case["window_history"] = history
    case["scheduler_metrics"] = metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--input-mode", choices=("static", "dynamic"), default="dynamic")
    parser.add_argument("--focused", action="store_true", help="run only the scheduler regressions")
    args = parser.parse_args()
    if args.focused and args.input_mode != "dynamic":
        parser.error("--focused requires --input-mode dynamic")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = build(work, args.input_mode)
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
        if args.input_mode == "dynamic":
            cases += [(2, 1, "queue", 2, 192, "decoded-refills", 64),
                      (3, 2, "queue", 2, 384, "window-reuse", 128),
                      (1, 4, "queue", 2, 24, "single-direct", 8),
                      (2, 2, "queue", 2, 16, "window-growth", 0),
                      (2, 4, "single", 2, 1, "pending-close", 0)]
        for ranks, ready, pipeline, mates, n, stress, member_records in cases:
            selected = records[:n]
            if stress == "window-growth":
                # With one record per chunk, only rank 0's third-window suffix grows.
                selected[10] = [(name, b"A" * 300000, b"I" * 300000)
                                for name, _, _ in selected[10]]
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
            output = tmp / "output.sam"
            env = dict(os.environ, SWBWA_GZIP_READY_CHUNKS=str(ready), HYDRA_LAUNCHER="fork")
            env.pop("SWBWA_OUTPUT_MD5", None)
            env.pop("HYDRA_HOST_FILE", None)
            env["PIPELINE_TEST_STRESS"] = stress
            command = ["mpirun", "-launcher", "fork", "-n", str(ranks), str(binary),
                       str(paths[0]), str(paths[1]) if mates == 2 else "-", str(output), pipeline]
            process = subprocess.run(command, env=env, text=True, capture_output=True, timeout=90)
            case = dict(ranks=ranks, ready=ready, pipeline=pipeline, mates=mates, records=n,
                        stress=stress, member_records=member_records,
                        returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
            logs.append(case)
            (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
            if stress == "pending-close":
                assert process.returncode != 0, case
                assert "MPI gzip close requires completed output batches" in process.stderr, case
                assert reports(process.stderr, "test pending close rank=0") == [
                    {"validated_EOF": "1", "pending": "1"}], case
                assert "pending close unexpectedly returned" not in process.stderr, case
                case["status"] = "EXPECTED_REJECTION"
                (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
                print("EXPECTED_REJECTION", stress, flush=True)
                continue
            assert process.returncode == 0, case
            assert output.read_bytes() == expected, case
            for rank in range(ranks):
                keys = reports(process.stderr, f"test pipeline keys rank={rank}")
                expected_pipelines = 2 + (pipeline == "queue")
                assert len(keys) == 1 and all(int(keys[0][key]) == expected_pipelines for key in (
                    "pipelines", "sets", "clears", "clear_direct_calls")), case
            if stress:
                check_stress(case)
            case["status"] = "PASS"
            (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
            output.unlink()
            print("PASS", ranks, ready, pipeline, mates, n, stress, flush=True)
    print("MPI_GZIP_ORDERED_PIPELINE_PASS", len(logs), flush=True)


if __name__ == "__main__":
    main()
