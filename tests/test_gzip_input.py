#!/usr/bin/env python3
"""Standalone host regression; does not edit/use the project Makefile."""
import argparse
import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "src" / "gzip"


def build(directory, sanitize, rapid=True):
    flags = ["-O1" if sanitize else "-O2", "-g", "-Wall", "-Wextra",
             "-D_FILE_OFFSET_BITS=64", "-pthread", "-I" + str(ROOT / "include"),
             "-I" + str(ENGINE), "-I" + str(ENGINE / "vendor" / "rapidgzip")]
    if rapid:
        flags += ["-DSWBWA_CODEC_RAPID=1", "-DSWBWA_CODEC_LDM=2"]
    if sanitize:
        flags += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer", "-fno-pie"]
    sources = [ENGINE / name for name in
               ("input.cpp", "gzip_stream.cpp", "rapid_reader.cpp", "codec_job.c", "backend_host.c")]
    if rapid:
        sources += [ENGINE / "rapid_firstpass.cpp"]
    sources += [ROOT / "tests" / "test_gzip_input.c"]
    objects = []
    for source in sources:
        obj = directory / (source.stem + ".o")
        cpp = source.suffix == ".cpp"
        compiler = os.environ.get("CXX" if cpp else "CC", "g++" if cpp else "gcc")
        standard = "c++17" if source.name == "rapid_firstpass.cpp" else "c++11" if cpp else "c11"
        subprocess.run([compiler, "-std=" + standard, *flags, "-c", str(source), "-o", str(obj)], check=True)
        objects.append(str(obj))
    executable = directory / "test-gzip-input"
    link = ["-fsanitize=address,undefined", "-no-pie"] if sanitize else []
    subprocess.run([os.environ.get("CXX", "g++"), *objects, *link,
                    "-lz", "-pthread", "-o", str(executable)], check=True)
    return executable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sanitize", action="store_true")
    parser.add_argument("--zlib-only", action="store_true")
    args = parser.parse_args()
    rng = random.Random(701)
    with tempfile.TemporaryDirectory(prefix="swbwa-gzip-input-") as temporary:
        directory = Path(temporary)
        executable = build(directory, args.sanitize, not args.zlib_only)
        env = dict(os.environ, SWBWA_GZIP_PREFETCH_MIB="0")
        passed = 0

        def case(name, first, second=None, *, encoded=None, encoded2=None,
                 error=False, chunk=65536, batch=None, workers=3, layout=0, backends=None):
            nonlocal passed
            plain = directory / "expected1.fq"
            plain2 = directory / "expected2.fq"
            source = directory / "input1.gz"
            source2 = directory / "input2.gz"
            plain.write_bytes(first)
            source.write_bytes(gzip.compress(first) if encoded is None else encoded)
            if second is not None:
                plain2.write_bytes(second)
                source2.write_bytes(gzip.compress(second) if encoded2 is None else encoded2)
            for backend in backends or (["zlib"] if args.zlib_only else ["rapid", "zlib"]):
                result = subprocess.run([str(executable), backend, str(source),
                    str(source2) if second is not None else "-", str(plain),
                    str(plain2) if second is not None else "-", str(int(error)), str(chunk),
                    str(workers), str(first.count(b"\n") // 4), str(layout), str(batch or chunk)],
                    env=env, text=True, capture_output=True, timeout=120)
                if result.returncode:
                    raise AssertionError(f"{name}/{backend}: {result.returncode}\n{result.stdout}\n{result.stderr}")
                print(f"PASS {name}/{backend}: {result.stdout.strip() or result.stderr.strip()}", flush=True)
                passed += 1

        def record(i, length=150, mate=1, comment=b""):
            sequence = bytes(rng.choice(b"ACGT") for _ in range(length))
            quality = bytes(rng.choice(b"!@+ABCDEFGHI") for _ in range(length))
            return b"@read%08d/%d" % (i, mate) + comment + b"\n" + sequence + b"\n+\n" + quality + b"\n"

        tiny = record(0, 5)
        records = [record(i, 40 + i % 191) for i in range(12000)]
        data = b"".join(records)
        same = data.replace(b"/1\n", b"/2\n")
        variable = b"".join(record(i, 10 + i % 53, 2, b" longer mate header") for i in range(12000))
        case("empty-se", b"")
        case("empty-pe", b"", b"")
        case("one-se", tiny)
        case("one-pe", tiny, tiny, layout=384)
        case("bulk-se", data)
        case("same-layout-pe", data, same, chunk=1 << 20, layout=384)
        case("same-layout-64cpe", data, same, chunk=1 << 20, layout=64, workers=1)
        case("independent-batch-size", data, same, chunk=262144, batch=2 << 20, layout=384)
        case("different-layout-reader", data, variable)
        case("different-layout-reversed", variable, data)
        case("different-batch-bytes-rejected", data, variable, chunk=1 << 20, layout=384, error=True)
        # Same total bytes AND equal records in the chunk, but shifted interior boundaries.
        fixed = b"".join(record(i, 50) for i in range(5000))
        shifted = record(0, 100) + b"".join(record(i, 50) for i in range(1, 4999)) + record(4999, 0)
        assert len(fixed) == len(shifted)
        case("equal-total-shifted-boundaries", fixed, shifted, chunk=1 << 20, layout=384, error=True)
        case("cross-chunk-long-record", record(0, 180000) + tiny, chunk=65536)
        # Four lines total exactly one reader chunk, including the final LF.
        exact = b"@xx\n" + b"A" * 32764 + b"\n+\n" + b"I" * 32764 + b"\n"
        assert len(exact) == 65536
        case("exact-chunk-eof", exact)
        case("missing-final-lf", tiny[:-1], error=True)
        case("missing-lf-exact-chunk", exact[:-1], error=True)
        for lines in (1, 2, 3):
            case(f"truncated-{lines}-lines", b"".join(tiny.splitlines(keepends=True)[:lines]), error=True)
        case("extra-partial-header", tiny + b"@truncated", error=True)
        case("one-empty-mate", tiny, b"", error=True)
        case("unequal-record-count", data, same + tiny, error=True)
        encoded = gzip.compress(data)
        case("truncated-gzip", data, encoded=encoded[:-4], error=True)
        corrupt = bytearray(encoded)
        corrupt[-8] ^= 1
        case("bad-crc", data, encoded=corrupt, error=True)
        corrupt = bytearray(encoded)
        corrupt[-4] ^= 1
        case("bad-isize", data, encoded=corrupt, error=True)
        case("raw-text-not-gzip", tiny, encoded=tiny, error=True)
        case("empty-file-not-gzip", b"", encoded=b"", error=True)
        case("trailing-garbage", tiny, encoded=gzip.compress(tiny) + b"garbage", error=True)
        case("zero-padding", tiny, encoded=gzip.compress(tiny) + bytes(1000))
        case("bad-padding", tiny, encoded=gzip.compress(tiny) + bytes(1000) + b"x", error=True)
        split = len(data) // 3 + 7
        chained = gzip.compress(data[:split]) + gzip.compress(b"") + gzip.compress(data[split:])
        case("concatenated-split-record", data, encoded=chained)
        case("stored-deflate", data, encoded=gzip.compress(data, compresslevel=0))
        compressor = zlib.compressobj(6, zlib.DEFLATED, 31, 8, zlib.Z_FIXED)
        case("fixed-deflate", data, encoded=compressor.compress(data) + compressor.flush())
        for size, workers in ((65537, 1), (131071, 6), (262144, 3)):
            case(f"odd-chunk-{size}-workers-{workers}", data, same, chunk=size, workers=workers)
        case("bad-chunk", tiny, chunk=1, error=True)
        case("bad-workers-zero", tiny, workers=0, error=True)
        case("bad-workers-over-limit", tiny, workers=385, error=True)
        case("bad-backend", tiny, backends=["auto"], error=True)
        if args.zlib_only:
            case("rapid-unavailable", tiny, backends=["rapid"], error=True)
        print(f"{passed} gzip input cases passed", flush=True)


if __name__ == "__main__":
    main()
