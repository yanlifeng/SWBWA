#!/usr/bin/env python3
"""Regression checks for the full ordered alignment comparator."""
from pathlib import Path
import tempfile
import unittest
from check_ordered_sam import compare, inspect


class OrderedSAMCheck(unittest.TestCase):
    def test_full_records_and_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq, reference, candidate = (root / x for x in ("reads.fq", "ref.sam", "out.sam"))
            fastq.write_bytes(b"@r1\nAC\n+\nII\n@r2\nGT\n+\nII\n")
            a = b"r1\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAC\tII\n"
            b = a.replace(b"r1", b"r2").replace(b"AC", b"GT")
            reference.write_bytes(b"@HD\tVN:1.6\n" + a + b)
            candidate.write_bytes(a + b)
            compare(reference, candidate)
            self.assertEqual(inspect(candidate, fastq, False)["reads"], 2)
            for wrong in (b + a, a, a + b + b, a + b.replace(b"GT", b"GG"),
                          a + b[:-1], a + b + b"@CO\tunexpected\n"):
                candidate.write_bytes(wrong)
                with self.assertRaises(ValueError):
                    compare(reference, candidate)
            for wrong in (b + a, a, a + b + b, a + b[:-1], a + b.replace(b"GT", b"G\0")):
                candidate.write_bytes(wrong)
                with self.assertRaises(ValueError):
                    inspect(candidate, fastq, False)

    def test_paired_supplementary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq, sam = root / "reads.fq", root / "out.sam"
            fastq.write_bytes(b"@r1/1\nAC\n+\nII\n")
            a = b"r1\t65\tchr1\t1\t60\t2M\t*\t0\t0\tAC\tII\n"
            b = a.replace(b"\t65\t", b"\t129\t")
            supplementary = a.replace(b"\t65\t", b"\t2113\t")
            sam.write_bytes(a + supplementary + b)
            self.assertEqual(inspect(sam, fastq, True)["records"], 3)
            sam.write_bytes(a + supplementary)
            with self.assertRaises(ValueError):
                inspect(sam, fastq, True)


if __name__ == "__main__":
    unittest.main()
