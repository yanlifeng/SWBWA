import copy
import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ldm_policy", ROOT / "tools/ldm_policy.py")
p = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = p
spec.loader.exec_module(p)


class PolicyTest(unittest.TestCase):
    def setUp(self):
        self.plan = {"training_dataset": "train79", "blocks": [["seed", "candidate", "seed", "seed"]] * 2,
                     "policies": {"seed": {}, "candidate": {}}}
        self.rows = [{"block": block, "order": order, "policy": name,
                      "dataset": "train79", "status": "PASS", "profile": 0,
                      "node": str(100 + block),
                      "digest": ["1", "20", "30", "40"],
                      "stage2": (100 if name == "seed" else 95) * (block + 1)}
                     for block, names in enumerate(self.plan["blocks"])
                     for order, name in enumerate(names)]

    def test_node_normalization(self):
        self.assertEqual(p.select(self.plan, self.rows)["selected"], "candidate")

    def test_completion_order_does_not_change_selection(self):
        self.assertEqual(p.select(self.plan, self.rows), p.select(self.plan, self.rows[::-1]))

    def test_noisy_incumbent_is_retained(self):
        self.rows[0]["stage2"] = 110
        self.assertEqual(p.select(self.plan, self.rows)["selected"], "seed")

    def test_incomplete_wrong_or_profiled_data_is_rejected(self):
        with self.assertRaises(ValueError):
            p.select(self.plan, self.rows[:-1])
        with self.assertRaises(ValueError):
            p.select(self.plan, self.rows + [self.rows[0]])
        for field, value in (("dataset", "validation"), ("status", "FAIL"),
                             ("profile", 1), ("stage2", float("nan")),
                             ("digest", ["wrong"])):
            rows = copy.deepcopy(self.rows)
            rows[0][field] = value
            with self.assertRaises(ValueError):
                p.select(self.plan, rows)
        rows = copy.deepcopy(self.rows)
        del rows[0]["profile"]
        with self.assertRaises(ValueError):
            p.select(self.plan, rows)

    def test_ownership_and_fixed_budget(self):
        for site in (0, 5, 8, 10):
            with self.assertRaises(ValueError):
                p.Policy().changed("caps", site, 256)
        for value in (-1, 32768, 1.5):
            with self.assertRaises(ValueError):
                p.Policy().changed("caps", 1, value)
        with self.assertRaises(ValueError):
            p.Policy(heap_cache=True)

    def test_size_statistics_do_not_choose_winner(self):
        counts = {name: [100, 10, 1, 0, 0, 0] for name in p.SITES}
        policies = p.candidates(counts, p.Policy(), "caps")
        self.assertIn("chain.caps.0", policies)
        self.assertIn("chain.caps.256", policies)
        self.assertEqual(p.quantile_cap([0, 0, 0, 0, 0, 1], .99), 16384)
        self.assertEqual(p.Policy().env()["SWBWA_LDM_PROFILE"], "0")

    def test_cold_start_has_no_ordinary_admission(self):
        seed = p.Policy(caps=(0,) * len(p.SITES))
        counts = {name: [100, 10, 1, 0, 0, 0] for name in p.SITES}
        policies = p.candidates(counts, seed, "admit")
        self.assertEqual(sum(policies["seed"].caps), 0)
        for name, value in policies.items():
            if name != "seed":
                self.assertEqual(sum(cap != 0 for cap in value.caps), 1)

    def test_mixed_node_block_is_rejected(self):
        self.rows[1]["node"] = "101"
        with self.assertRaises(ValueError):
            p.select(self.plan, self.rows)


if __name__ == "__main__":
    unittest.main()
