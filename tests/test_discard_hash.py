import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "check_discard_hash", Path(__file__).resolve().parents[1] /
    "scripts/check_discard_hash.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class DiscardHashTests(unittest.TestCase):
    def parse(self, lines):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.log"
            path.write_text(lines)
            return checker.parse_log(path)

    def line(self, rank, suffix=""):
        return (f"[SWBWA output hash rank {rank:06d}/000002] "
                "calls=1 bytes=20 sum=0x0000000000000003 "
                f"xor=0x0000000000000003 enabled=1{suffix}\n")

    def test_legacy_and_full_hash(self):
        for suffix in ("", " hash_prefix_bytes=0"):
            self.assertEqual(self.parse(self.line(0, suffix) + self.line(1, suffix)),
                             (2, 2, 40, 6, 0))

    def test_reject_prefix_hash(self):
        with self.assertRaisesRegex(ValueError, "prefix-only"):
            self.parse(self.line(0, " hash_prefix_bytes=64") + self.line(1))

    def test_reject_missing_duplicate_and_disabled(self):
        for lines in (self.line(0), self.line(0) * 2,
                      (self.line(0) + self.line(1)).replace("enabled=1", "enabled=0")):
            with self.assertRaises(ValueError):
                self.parse(lines)


if __name__ == "__main__":
    unittest.main()
