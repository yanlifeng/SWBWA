"""Scalar algebra model for fused CPE gap updates, not a Sunway execution test.

The production helpers use signed int32 differences and 0/1 comparison lanes.
This checks their select/clamp algebra against the original separate clamps;
target compilation and whole-alignment fingerprints are still required.
"""

import itertools
import random
import unittest


def original(gap, extension, h, open_extension):
    return max(max(gap - extension, 0), max(h - open_extension, 0))


def fused(gap, extension, h, open_extension):
    a, b = gap - extension, h - open_extension
    mask = -int(a < b)
    value = (mask & b) | ((mask ^ -1) & a)
    return value & -int(0 < value)


def fused_lanes(gap, extension, h, open_extension, lanes):
    return [fused(*values) if i < lanes else 0
            for i, values in enumerate(zip(gap, extension, h, open_extension))]


class FusedGapTests(unittest.TestCase):
    def test_u8_exhaustive_difference_pairs(self):
        # Every difference of two unsigned bytes, including zero/equal ties.
        for a, b in itertools.product(range(-255, 256), repeat=2):
            args = (max(a, 0), max(-a, 0), max(b, 0), max(-b, 0))
            self.assertEqual(fused(*args), original(*args), args)

    def test_i16_boundaries(self):
        values = (-32768, -32767, -255, -1, 0, 1, 255, 32766, 32767)
        for args in itertools.product(values, repeat=4):
            self.assertEqual(fused(*args), original(*args), args)

    def test_random_mixed_lanes(self):
        rng = random.Random(7963242)
        for lanes in (8, 16):
            for case in range(4096):
                # Also exercise sign-extended penalty broadcasts and dirty
                # inactive i16 lanes; only the low eight may survive.
                lo, hi = (-32768, 32767) if lanes == 8 else (-128, 255)
                vectors = [[rng.randint(lo, hi) for _ in range(16)]
                           for _ in range(4)]
                expected = [original(*args) if i < lanes else 0
                            for i, args in enumerate(zip(*vectors))]
                self.assertEqual(fused_lanes(*vectors, lanes), expected,
                                 (lanes, case))

    def test_dependent_gap_updates(self):
        rng = random.Random(16)
        for limit in (255, 32767):
            for _ in range(100):
                old = new = 0
                extension, opening = rng.randint(0, 20), rng.randint(0, 40)
                for _ in range(200):
                    h = rng.randint(0, limit)
                    old = original(old, extension, h, opening)
                    new = fused(new, extension, h, opening)
                    self.assertEqual(new, old)


if __name__ == "__main__":
    unittest.main()
