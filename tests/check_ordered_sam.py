#!/usr/bin/env python3
"""Check complete SAM contents and FASTQ order, without sorting or rewriting SAM."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def fastq_names(path):
    with Path(path).open("rb") as source:
        while True:
            header = source.readline()
            if not header:
                return
            seq, plus, qual = (source.readline() for _ in range(3))
            if (not header.startswith(b"@") or not plus.startswith(b"+") or
                    not qual or len(seq.rstrip(b"\r\n")) != len(qual.rstrip(b"\r\n"))):
                raise ValueError("Invalid four-line FASTQ")
            name = header[1:].split()[0]
            if name.endswith((b"/1", b"/2")):
                name = name[:-2]
            yield name


def sam_lines(path):
    in_header = True
    prefixes = (b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t")
    with Path(path).open("rb") as source:
        for line in source:
            if line.startswith(b"@"):
                if not in_header or not line.startswith(prefixes):
                    raise ValueError("Invalid or misplaced SAM header in %s" % path)
                continue
            in_header = False
            yield line


def inspect(path, fastq, paired):
    names = iter(fastq_names(fastq))
    expected = next(names, None)
    primary = []
    records = reads = size = 0
    digest = hashlib.md5()

    def finish_read():
        required = [64, 128] if paired else [0]
        if sorted(primary) != required:
            raise ValueError("Missing/duplicate primary alignment for %r: %r" %
                             (expected, primary))

    for line in sam_lines(path):
        fields = line.rstrip(b"\n").split(b"\t")
        if b"\0" in line or not line.endswith(b"\n") or len(fields) < 11:
            raise ValueError("Malformed SAM record %d in %s" % (records + 1, path))
        if fields[0] != expected:
            finish_read()
            reads += 1
            expected = next(names, None)
            primary = []
        if fields[0] != expected:
            raise ValueError("SAM is not in FASTQ order: expected %r, got %r" %
                             (expected, fields[0]))
        flag = int(fields[1])
        if not flag & (0x100 | 0x800):
            if bool(flag & 1) != paired:
                raise ValueError("Incorrect paired flag")
            primary.append(flag & (64 | 128))
        digest.update(line)
        size += len(line)
        records += 1
    if expected is not None:
        finish_read()
        reads += 1
        if next(names, None) is not None:
            raise ValueError("SAM is missing trailing FASTQ records")
    return dict(records=records, reads=reads, bytes=size, md5=digest.hexdigest())


def compare(reference, candidate):
    for number, (left, right) in enumerate(itertools.zip_longest(
            sam_lines(reference), sam_lines(candidate)), 1):
        if left != right:
            raise ValueError("Full SAM mismatch at record %d: %s versus %s" %
                             (number, reference, candidate))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--candidate", action="append", required=True)
    args = parser.parse_args()
    reference = inspect(args.reference, args.fastq, args.paired)
    print("REFERENCE " + json.dumps(dict(path=args.reference, **reference)), flush=True)
    for path in args.candidate:
        compare(args.reference, path)
        actual = inspect(path, args.fastq, args.paired)
        if actual != reference:
            raise ValueError("SAM statistics differ")
        print("ORDERED_SAM_PASS " + json.dumps(dict(path=path, **actual)), flush=True)


if __name__ == "__main__":
    main()
