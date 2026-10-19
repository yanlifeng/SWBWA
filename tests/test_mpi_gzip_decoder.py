#!/usr/bin/env python3
"""Standalone local host-backend MPI tests; never invokes a remote launcher."""
import argparse
import gzip
import json
import os
from pathlib import Path
import random
import shlex
import struct
import subprocess
import tempfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "src" / "gzip"


def build(work, sanitizer):
    flags = ["-O1" if sanitizer else "-O2", "-g", "-Wall", "-Wextra",
             "-D_FILE_OFFSET_BITS=64", "-D_POSIX_C_SOURCE=200809L", "-DSWBWA_CODEC_RAPID=1", "-DSWBWA_CODEC_LDM=2",
             "-I" + str(ROOT / "include"), "-I" + str(ENGINE), "-I" + str(ENGINE / "vendor" / "rapidgzip")]
    if sanitizer:
        flags += ["-fsanitize=" + sanitizer, "-fno-sanitize-recover=all", "-fno-omit-frame-pointer"]
        if "address" in sanitizer:
            flags += ["-fno-pie"]
    sources = [ENGINE / n for n in ("mpi_decoder.cpp", "codec_job.c", "backend_host.c", "rapid_firstpass.cpp")]
    sources += [ROOT / "tests" / "test_mpi_gzip_decoder.cpp"]
    objects = []
    for source in sources:
        obj = work / (source.stem + ".o")
        cpp = source.suffix == ".cpp"
        compiler = os.environ.get("MPICXX", "mpicxx") if cpp else os.environ.get("CC", "cc")
        standard = "c++17" if source.name == "rapid_firstpass.cpp" else "c++11" if cpp else "c11"
        subprocess.run([compiler, "-std=" + standard, *flags, "-c", str(source), "-o", str(obj)], check=True)
        objects.append(str(obj))
    binary = work / "test-mpi-gzip-decoder"
    subprocess.run([os.environ.get("MPICXX", "mpicxx"), *objects,
                    *(["-fsanitize=" + sanitizer] if sanitizer else []),
                    *(["-no-pie"] if sanitizer and "address" in sanitizer else []),
                    "-Wl,--wrap=swbwa_codec_launch", "-Wl,--wrap=MPI_Comm_free", "-Wl,--wrap=fstat64",
                    "-Wl,--wrap=__fxstat64",
                    "-lz", "-pthread", "-o", str(binary)], check=True)
    return binary


def wrap(raw, body, header=b"\x1f\x8b\x08\0\0\0\0\0\0\x03"):
    return header + body + struct.pack("<II", zlib.crc32(raw), len(raw) & 0xffffffff)


def encoded(raw, strategy=zlib.Z_DEFAULT_STRATEGY, level=6):
    codec = zlib.compressobj(level, zlib.DEFLATED, -15, 8, strategy)
    return wrap(raw, codec.compress(raw) + codec.flush())


def stored(raw, block=65531):
    parts = []
    for i in range(0, len(raw), block):
        part = raw[i:i + block]
        parts.append(b"\0" + struct.pack("<HH", len(part), len(part) ^ 65535) + part)
    return wrap(raw, b"".join(parts) + b"\x01\0\0\xff\xff")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--sanitize", choices=("undefined", "address,undefined"))
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve() if args.binary else build(work, args.sanitize)
    subprocess.run([str(binary), "--unit"], check=True)
    # fork is explicitly local, independent of site MPI launcher defaults.
    launcher = ["mpirun", "-launcher", "fork"]
    logs, reports = [], []
    rng = random.Random(7315)
    with tempfile.TemporaryDirectory(prefix="fixtures-", dir=work) as temporary:
        root = Path(temporary)

        def run(name, raw, blob=None, ranks=3, jobs=4, capacity=1 << 20,
                chunk=32768, lookahead=1 << 20, memory=32 << 20,
                decoded=8 << 20, mode="normal", error=None, workers=64):
            source, reference = root / (name + ".gz"), root / (name + ".raw")
            source.write_bytes(encoded(raw) if blob is None else blob)
            reference.write_bytes(raw)
            command = launcher + ["-n", str(ranks), str(binary), str(source), str(reference),
                *map(str, (jobs, capacity, chunk, lookahead, memory, decoded)), mode, str(workers)]
            process = subprocess.run(command, text=True, capture_output=True, timeout=90)
            logs.append(dict(case=name, command=shlex.join(command), returncode=process.returncode,
                             stdout=process.stdout, stderr=process.stderr))
            if error:
                assert process.returncode and error in process.stderr, logs[-1]
                assert "MPI_GZIP_RESULT" not in process.stdout, logs[-1]
                reports.append(dict(case=name, status="EXPECTED_REJECTION"))
                print("EXPECTED_REJECTION", name, flush=True)
                return None
            assert process.returncode == 0, logs[-1]
            lines = [s.removeprefix("MPI_GZIP_RESULT ") for s in process.stdout.splitlines()
                     if s.startswith("MPI_GZIP_RESULT ")]
            assert len(lines) == 1, logs[-1]
            result = json.loads(lines[0])
            assert result["decoded_bytes"] == len(raw) and result["crc32"] == zlib.crc32(raw)
            assert sum(s[4] for s in result["stats"]) == len(raw)
            assert all(s[6] <= memory for s in result["stats"])
            reports.append(dict(case=name, status="PASS", result=result))
            print("PASS", name, "epochs", result["epochs"], flush=True)
            return result

        try:
            run("empty", b"")
            run("tiny-more-ranks-than-jobs", b"@a\nAC\n+\nII\n", ranks=5, jobs=1)
            run("close-mpi-error", b"close\n", mode="close-error", error="decoder close failed:")
            dna = bytes(rng.choices(b"ACGTN\n", weights=(25, 25, 25, 25, 1, 2), k=2 << 20))
            dynamic = encoded(dna)
            for ranks in (1, 2, 3, 5):
                result = run(f"dynamic-ranks{ranks}", dna, dynamic, ranks=ranks)
                assert max(s[5] for s in result["stats"]) > 1
                assert sum(s[1] for s in result["stats"]) > 0
                assert all(s[0] > 0 and s[4] > 0 for s in result["stats"])
            run("reversed-communicator", dna, dynamic, mode="reverse")
            run("independent-dual-readers", dna, dynamic, mode="dual")
            for jobs in (32, 64, 384):
                run(f"jobs{jobs}", dna, dynamic, jobs=jobs, capacity=131072, chunk=4096, workers=384)
            run("default-budgets", dna, dynamic, jobs=64, capacity=2 << 20, chunk=256 << 10,
                memory=512 << 20, decoded=128 << 20, workers=384)
            random_bytes = rng.randbytes((1 << 20) + 17)
            run("stored", random_bytes, stored(random_bytes))
            result = run("full-384-job-launch", random_bytes, stored(random_bytes), ranks=1,
                jobs=384, capacity=131072, chunk=1024, memory=256 << 20, decoded=64 << 20, workers=384)
            assert result["stats"][0][5] == 384
            run("short-history", random_bytes[:10000], stored(random_bytes[:10000], 37),
                chunk=1024, lookahead=16384)
            planted = bytearray(random_bytes)
            planted[32763:32777] = b"\0\x03\0\xfc\xffABC\x01\x01\0\xfe\xffZ"
            result = run("false-candidate", bytes(planted), stored(bytes(planted)))
            assert sum(s[2] for s in result["stats"]) > 0
            assert sum(w["discarded_output_bytes"] for w in result["work"]) > 0
            engine = zlib.compressobj(6, zlib.DEFLATED, -15)
            parts, blocks, previous = [], [], dna[:32768]
            for i in range(24):
                part = previous[-32768:] * 2 + dna[i * 32768:(i + 1) * 32768]
                parts.append(part); blocks.append(engine.compress(part) + engine.flush(zlib.Z_BLOCK))
                previous = part
            blocks.append(engine.flush())
            dependency = b"".join(parts)
            run("cross-job-history", dependency, wrap(dependency, b"".join(blocks)))
            header = bytearray(b"\x1f\x8b\x08\x1e\0\0\0\0\0\x03")
            header += b"\x04\0abcdreads.fastq\0comment\0"
            header += struct.pack("<H", zlib.crc32(header) & 65535)
            optional = bytes(header) + dynamic[10:]
            run("optional-header", dna, optional)
            for level in (0, 1, 6, 9):
                run(f"gzip-level{level}", dna, encoded(dna, level=level))
            for name, strategy in (("fixed", zlib.Z_FIXED), ("huffman", zlib.Z_HUFFMAN_ONLY), ("rle", zlib.Z_RLE)):
                run(name, dna, encoded(dna, strategy), lookahead=4 << 20)
            high = b"A\n" * (2 << 20)
            result = run("capacity-retry", high, encoded(high, zlib.Z_FIXED), capacity=131072,
                         decoded=8 << 20, memory=32 << 20, jobs=64)
            assert all(s[3] >= 1 for s in result["stats"])
            assert len({s[3] for s in result["stats"]}) == 1
            assert sum(w["retry_jobs"] for w in result["work"]) > 0
            assert sum(w["discover_output_full_jobs"] for w in result["work"]) > 0
            assert all(t["retry_s"] > 0 for t in result["timers"])
            irregular = b"T\n" * (1280 << 10)
            run("non-power-of-two-budget", irregular, encoded(irregular, zlib.Z_FIXED), capacity=131072,
                decoded=3 << 20, memory=12 << 20, jobs=64)
            run("odd-chunk", dna, dynamic, chunk=65537)
            run("one-worker-many-jobs", dna, dynamic, workers=1)
            run("members", dna + b"extra", dynamic + encoded(b"extra"))
            run("empty-members", b"", encoded(b"") * 5, jobs=2)
            chained = encoded(b"") + encoded(dna[:100003]) + encoded(b"") + encoded(dna[100003:]) + encoded(b"")
            run("mixed-empty-members", dna, chained, jobs=2, chunk=16384)
            run("members-reset-dictionary", dna + random_bytes, dynamic + stored(random_bytes), jobs=2)
            run("timestamp-change-during-decode", dna, dynamic, jobs=1, mode="touch")
            run("rank-timestamp-skew", dna, dynamic, mode="stat-time-skew")
            run("config-mismatch", dna, dynamic, mode="config-mismatch",
                error="metadata/configuration: chunk_bytes local=")
            run("path-mismatch", dna, dynamic, mode="path-mismatch", error="input path")
            run("invalid-one-rank", dna, dynamic, mode="invalid-local", error="invalid decoder configuration")
            corrupt = bytearray(dynamic); corrupt[-8] ^= 1
            run("bad-crc", dna, corrupt, error="CRC/ISIZE mismatch")
            corrupt = bytearray(dynamic); corrupt[-4] ^= 1
            run("bad-isize", dna, corrupt, error="CRC/ISIZE mismatch")
            corrupt = bytearray(optional); corrupt[len(header) - 1] ^= 1
            run("bad-header-crc", dna, corrupt, error="header CRC mismatch")
            corrupt = bytearray(encoded(b"extra")); corrupt[-8] ^= 1
            run("bad-second-member-crc", dna, dynamic + corrupt, error="CRC/ISIZE mismatch")
            run("bad-second-member-header", dna, dynamic + b"not-a-gzip-member", error="invalid gzip header")
            prefix = random_bytes[:32768]
            dependent = prefix[-16384:]
            engine = zlib.compressobj(9, zlib.DEFLATED, -15, zdict=prefix)
            illegal = wrap(dependent, engine.compress(dependent) + engine.flush())
            run("illegal-cross-member-dictionary", prefix + dependent, encoded(prefix) + illegal,
                error="MPI_GZIP_ERROR")
            run("padding", dna, dynamic + bytes(16), error="invalid gzip header")
            run("truncated-trailer", dna, dynamic[:-4], error="truncated trailer")
            run("truncated-body", dna, dynamic[:len(dynamic) // 2], error="anchored DISCOVER failed")
            run("truncated-header", b"", b"\x1f\x8b\x08\x08" + bytes(6) + b"x" * 30,
                error="truncated gzip header")
            run("reserved-flags", dna, dynamic[:3] + b"\x20" + dynamic[4:], error="reserved gzip")
            corrupt = bytearray(dynamic); corrupt[10] = (corrupt[10] & ~7) | 6
            run("invalid-deflate", dna, corrupt, error="anchored DISCOVER failed")
            corrupt = bytearray(stored(random_bytes)); corrupt[13] ^= 1
            run("invalid-stored-length", random_bytes, corrupt, error="anchored DISCOVER failed")
            run("empty-file", b"", b"", error="regular gzip file")
            run("output-budget", high, encoded(high, zlib.Z_FIXED), capacity=131072,
                decoded=1 << 20, error="exhausted memory/decoded budget")
            run("input-lookahead-budget", random_bytes, stored(random_bytes), chunk=1024, lookahead=1024,
                error="insufficient lookahead")
        finally:
            (work / "test-results.json").write_text(json.dumps(reports, indent=2) + "\n")
            (work / "test-logs.json").write_text(json.dumps(logs, indent=2) + "\n")
    print("MPI_GZIP_DECODER_TESTS_PASS", len(reports), flush=True)


if __name__ == "__main__":
    main()
