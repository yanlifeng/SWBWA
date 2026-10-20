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


def build(work, sanitizer, input_mode="static"):
    flags = ["-O1" if sanitizer else "-O2", "-g", "-Wall", "-Wextra",
             "-D_FILE_OFFSET_BITS=64", "-D_POSIX_C_SOURCE=200809L",
             "-DSWBWA_CODEC_RAPID=1", "-DSWBWA_CODEC_LDM=2", "-DSWBWA_USE_MPI=1",
             "-DSWBWA_MPI_INPUT_MODE=SWBWA_MPI_INPUT_" + input_mode.upper(),
             "-I" + str(ROOT / "include"), "-I" + str(ENGINE),
             "-I" + str(ENGINE / "vendor" / "rapidgzip")]
    san = ["-fsanitize=" + sanitizer, "-fno-sanitize-recover=all", "-fno-omit-frame-pointer"] if sanitizer else []
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
    subprocess.run([os.environ.get("MPICXX", "mpicxx"), *objects, *san,
                    "-Wl,--wrap=swbwa_mpi_decoder_config_init", "-Wl,--wrap=swbwa_mpi_decoder_next",
                    "-lz", "-pthread", "-o", str(binary)], check=True)
    frame = work / "test-fastq-frame"
    subprocess.run([os.environ.get("CXX", "c++"), "-std=c++11", *flags,
                    str(ROOT / "tests" / "test_fastq_frame.cpp"), str(work / "fastq_frame.o"),
                    "-o", str(frame)], check=True)
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
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.binary:
        binary, frame = args.binary.resolve(), work / "test-fastq-frame"
    else:
        binary, frame = build(work, args.sanitize, args.input_mode)
    subprocess.run([str(frame)], check=True)
    rng = random.Random(2817)
    reports, logs = [], []
    with tempfile.TemporaryDirectory(prefix="fixtures-", dir=work) as temporary:
        root = Path(temporary)

        def run(name, r1, r2=None, blob1=None, blob2=None, ranks=2,
                target=32768, mode="normal", error=None, defaults=False, lookahead=1 << 20, ready=4):
            raw1, gz1 = root / (name + ".R1.raw"), root / (name + ".R1.gz")
            raw1.write_bytes(r1); gz1.write_bytes(encoded(r1) if blob1 is None else blob1)
            raw2 = gz2 = "-"
            if r2 is not None:
                raw2, gz2 = root / (name + ".R2.raw"), root / (name + ".R2.gz")
                raw2.write_bytes(r2); gz2.write_bytes(encoded(r2) if blob2 is None else blob2)
            command = ["mpirun", "-launcher", "fork", "-n", str(ranks), str(binary),
                       str(gz1), str(gz2), str(raw1), str(raw2), str(target), "384", "1", mode]
            env = dict(os.environ)
            env["FASTQ_TEST_LOOKAHEAD"] = str(lookahead)
            env["SWBWA_GZIP_READY_CHUNKS"] = str(ready)
            if defaults:
                env["FASTQ_TEST_DEFAULT_CONFIG"] = "1"
            else:
                env.pop("FASTQ_TEST_DEFAULT_CONFIG", None)
            env.setdefault("ASAN_OPTIONS", "detect_leaks=0")
            process = subprocess.run(command, text=True, capture_output=True, timeout=120, env=env)
            logs.append(dict(case=name, command=shlex.join(command), returncode=process.returncode,
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
            estimates = set(map(int, re.findall(r"records_per_chunk=(\d+)", process.stderr)))
            assert len(estimates) == 1, logs[-1]
            k = estimates.pop()
            for rank, chunk, first, end, count in re.findall(
                r"\[MPI gzip rank=(\d+)\] chunk=(\d+) first_record=(\d+) end_record_exclusive=(\d+) records=(\d+)", process.stderr
            ):
                rank, chunk, first, end, count = map(int, (rank, chunk, first, end, count))
                assert end == first + count
                if args.input_mode == "static":
                    assert chunk % ranks == rank
                assert first == min(chunk * k, result["records"])
            result["boundary_records"] = sum(map(int, re.findall(r"boundary_records=(\d+)", process.stderr)))
            result["max_epochs"] = max(map(int, re.findall(r"epochs=(\d+)", process.stderr)))
            result["consumer_chunks"] = dict((int(rank), int(count)) for rank, count in re.findall(
                r"\[test consumer rank=(\d+)\] chunks=(\d+)", process.stderr))
            result["remote_claims"] = sum(map(int, re.findall(r"remote_claims=(\d+)", process.stderr)))
            reports.append(dict(case=name, status="PASS", result=result))
            print("PASS", name, result, flush=True)
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
            run("default-decoder-config", r1, r2, defaults=True)
            run("thread-multiple-reader", r1, r2, mode="thread")
            run("single-end", r1, ranks=3)
            for name, strategy in (("fixed", zlib.Z_FIXED), ("huffman", zlib.Z_HUFFMAN_ONLY), ("rle", zlib.Z_RLE)):
                run(name, r1, r2, blob1=encoded(r1, strategy), blob2=encoded(r2, strategy), lookahead=4 << 20)
            run("stored", r1, r2, blob1=stored(r1), blob2=stored(r2))
            concat1 = encoded(b"") + encoded(r1[:43]) + encoded(r1[43:451231]) + encoded(b"") + encoded(r1[451231:])
            concat2 = encoded(r2[:10101]) + encoded(b"") + encoded(r2[10101:]) + encoded(b"")
            run("concat-members-cross-record", r1, r2, blob1=concat1, blob2=concat2, ranks=3)
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
                for ready in (1, 2, 8):
                    run(f"ready-window-{ready}", r1, r2, ranks=3, ready=ready)
                result = run("skewed-consumer", r1, r2, ranks=3, mode="skew", ready=8)
                assert result["remote_claims"] > 0
                counts = result["consumer_chunks"]
                assert counts[0] + counts[1] > 2 * counts[2], counts
                run("bad-ready-limit", tiny, ready=0, mode="reject", error="SWBWA_GZIP_READY_CHUNKS")
                run("ready-config-mismatch", tiny, mode="window-mismatch", error="reader configuration mismatch")
                run("invalid-completion", tiny, mode="bad-complete", error="completion out of order")
        finally:
            (work / "test-results.json").write_text(json.dumps(reports, indent=2) + "\n")
            (work / "test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
    print("MPI_FASTQ_INPUT_TESTS_PASS", len(reports), flush=True)


if __name__ == "__main__":
    main()
