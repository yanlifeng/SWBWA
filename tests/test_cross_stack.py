#!/usr/bin/env python3
"""Validate -v 4 cross-runtime logs, or run parser tests with no arguments.

Hardware regression: use the same input and -K for cgs and cgs_cross builds
with OUTPUT_MODE=discard and DISCARD_HASH_BYTES=0, then compare the complete
output hash (calls/bytes/sum/xor). Choose -K to produce at least 300 batches.
This checker validates dispatch/stack traces, not SAM correctness by itself.
"""
import re
import sys
import unittest
from pathlib import Path

TRACE = re.compile(r"\[cross-stack\] task=(\d+) sp_min=(0x[0-9a-f]+) sp_max=(0x[0-9a-f]+)")


def check_trace(text, min_batches=1):
    rows = [(int(n), int(lo, 16), int(hi, 16)) for n, lo, hi in TRACE.findall(text)]
    if not rows:
        raise ValueError("No cross-stack trace; run cgs_cross with -v 4")
    expected = rows[0][1:]
    if not 0 < expected[0] <= expected[1] <= 256 * 1024 or expected[0] % 16:
        raise ValueError("Invalid stack range")
    if expected[0] != expected[1]:
        raise ValueError("Nonuniform CPE stacks; extrema alone cannot prove per-core stability")
    for index, (sequence, lo, hi) in enumerate(rows, 1):
        if sequence != index:
            raise ValueError(f"Missing or duplicate task at {index}: {sequence}")
        if (lo, hi) != expected:
            raise ValueError(f"Stack drift at task {sequence}: {(lo, hi)} != {expected}")
    batches = len(re.findall(r"\[M::mem_process_seqs_merge2\] Processed \d+ reads", text))
    if batches == 0 or len(rows) != 3 * batches:
        raise ValueError("Incomplete three-phase batches")
    if batches < min_batches:
        raise ValueError(f"Only {batches} batches; require at least {min_batches}")
    if "[main] CMD:" not in text:
        raise ValueError("Missing normal program completion")
    if "UNIT_SCORE exception" in text or "[E::" in text:
        raise ValueError("Run contains an error")
    return batches, len(rows), expected


class CrossStackTest(unittest.TestCase):
    @staticmethod
    def log():
        return "\n".join([
            f"[cross-stack] task={i} sp_min=0x3ff80 sp_max=0x3ff80" for i in range(1, 7)
        ] + ["[M::mem_process_seqs_merge2] Processed 12 reads"] * 2 + ["[main] CMD: SWBWA mem"])

    def test_stable(self):
        self.assertEqual(check_trace(self.log()), (2, 6, (0x3ff80, 0x3ff80)))

    def test_drift(self):
        with self.assertRaisesRegex(ValueError, "drift"):
            check_trace(self.log().replace("task=6 sp_min=0x3ff80", "task=6 sp_min=0x3fd80"))

    def test_missing_dispatch(self):
        with self.assertRaisesRegex(ValueError, "Missing or duplicate"):
            check_trace(self.log().replace("task=4", "task=5"))

    def test_incomplete(self):
        with self.assertRaisesRegex(ValueError, "completion"):
            check_trace(self.log().replace("[main] CMD:", "interrupted:"))

    def test_invalid_stack(self):
        with self.assertRaisesRegex(ValueError, "Invalid stack"):
            check_trace(self.log().replace("0x3ff80", "0x1"))

    def test_nonuniform_stack(self):
        with self.assertRaisesRegex(ValueError, "Nonuniform"):
            check_trace(self.log().replace("sp_max=0x3ff80", "sp_max=0x3ffc0"))

    def test_stress_threshold(self):
        with self.assertRaisesRegex(ValueError, "require at least 300"):
            check_trace(self.log(), min_batches=300)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        unittest.main()
    else:
        import argparse
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--min-batches", type=int, default=1)
        parser.add_argument("logs", nargs="+")
        args = parser.parse_args()
        for name in args.logs:
            batches, tasks, (low, high) = check_trace(
                Path(name).read_text(errors="replace"), args.min_batches)
            print(f"PASS {name}: {batches} batches, {tasks} tasks, SP=[{low:#x}, {high:#x}]")
