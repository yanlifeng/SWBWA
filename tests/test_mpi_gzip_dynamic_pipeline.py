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
                    "-pthread", "-o", str(binary)], check=True)
    return binary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--input-mode", choices=("static", "dynamic"), default="dynamic")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = build(work, args.input_mode)
    rng = random.Random(509)
    records = []
    for i in range(227):
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
        cases = [(ranks, ready, pipeline, mates, n)
                 for ranks in (1, 2, 3, 6) for pipeline in ("single", "queue")
                 for mates in (1, 2) for ready in (1, 4) for n in (227,)]
        cases += [(3, 4, pipeline, 2, n) for pipeline in ("single", "queue") for n in (0, 1)]
        for ranks, ready, pipeline, mates, n in cases:
            paths = [tmp / "R1.gz", tmp / "R2.gz"]
            for m, path in enumerate(paths):
                raw = b"".join(b"@" + pair[m][0] + b" comment\n" + pair[m][1] +
                               b"\n+\n" + pair[m][2] + b"\n" for pair in records[:n])
                path.write_bytes(gzip.compress(raw, mtime=0))
            expected = b"".join(pair[m][0] + b"\t4\t*\t0\t0\t*\t*\t0\t0\t" + pair[m][1] +
                                b"\t" + pair[m][2] + f"\tGI:i:{i * mates + m}\n".encode()
                                for i, pair in enumerate(records[:n]) for m in range(mates))
            output = tmp / "output.sam"
            env = dict(os.environ, SWBWA_GZIP_READY_CHUNKS=str(ready), HYDRA_LAUNCHER="fork")
            env.pop("SWBWA_OUTPUT_MD5", None)
            env.pop("HYDRA_HOST_FILE", None)
            command = ["mpirun", "-launcher", "fork", "-n", str(ranks), str(binary),
                       str(paths[0]), str(paths[1]) if mates == 2 else "-", str(output), pipeline]
            process = subprocess.run(command, env=env, text=True, capture_output=True, timeout=90)
            case = dict(ranks=ranks, ready=ready, pipeline=pipeline, mates=mates, records=n,
                        returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
            logs.append(case)
            (work / "pipeline-test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
            assert process.returncode == 0, case
            assert output.read_bytes() == expected, case
            output.unlink()
            print("PASS", ranks, ready, pipeline, mates, n, flush=True)
    print("MPI_GZIP_ORDERED_PIPELINE_PASS", len(logs), flush=True)


if __name__ == "__main__":
    main()
