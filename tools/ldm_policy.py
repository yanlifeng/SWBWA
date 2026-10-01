#!/usr/bin/env python3
"""Bounded, offline admission-policy search; never migrate live allocations.

Profiles propose sizes, not access heat. Only measured, correctness-checked
Stage2 improvements can change the incumbent. Ordinary admission is tuned;
ownership rules, manual scratch hints and allocator representation stay fixed.
"""
import argparse
import dataclasses
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics

SITES = ("other", "chain", "reference", "global_dp", "reg2aln", "query_dp",
         "extend_dp", "smem", "chain_seed", "context", "dedup_sort")
TUNABLE = (1, 2, 3, 4)
BOUNDS = (256, 1024, 4096, 8192, 16384)
PAYLOAD = 56 * 1024


@dataclasses.dataclass(frozen=True)
class Policy:
    caps: tuple = (0, 4096, 1024, 4096, 256, 0, 0, 0, 0, 0, 0)
    reserves: tuple = (0, 4096, 4096, 0, 4096, 0, 0, 0, 0, 0, 0)
    hints: int = 255
    heap_cache: int = 1
    pool_bitmap: int = 2

    def __post_init__(self):
        if len(self.caps) != len(SITES) or len(self.reserves) != len(SITES):
            raise ValueError("Wrong site schema")
        for values in (self.caps, self.reserves):
            if any(type(x) is not int or x < 0 or x > PAYLOAD for x in values):
                raise ValueError("Invalid byte limit")
        if self.caps[0] or any(self.caps[i] for i in range(5, len(SITES))):
            raise ValueError("Search cannot admit unaudited or hint-managed sites")
        if any(self.caps[i] > 16384 for i in TUNABLE):
            raise ValueError("Search is bounded to small ordinary objects")
        fixed = (self.hints, self.heap_cache, self.pool_bitmap)
        if any(type(x) is not int for x in fixed) or fixed != (255, 1, 2):
            raise ValueError("Search must keep scratch/metadata configuration fixed")

    def changed(self, field, site, value):
        values = list(getattr(self, field))
        values[site] = value
        return dataclasses.replace(self, **{field: tuple(values)})

    def env(self):
        return {"SWBWA_LDM_CAPS": ",".join(map(str, self.caps)),
                "SWBWA_LDM_RESERVES": ",".join(map(str, self.reserves)),
                "SWBWA_LDM_SCRATCH_HINTS": str(self.hints),
                "SWBWA_LDM_POOL_CACHE": str(self.heap_cache),
                "SWBWA_LDM_POOL_BITMAP": str(self.pool_bitmap),
                "SWBWA_LDM_GROW_IN_PLACE": "1", "SWBWA_LDM_PROFILE": "0"}


def from_dict(value):
    value = dict(value)
    value["caps"] = tuple(value["caps"])
    value["reserves"] = tuple(value["reserves"])
    return Policy(**value)


def histogram(text):
    result = {name: [0] * 6 for name in SITES}
    seen = set()
    pattern = r"\[LDM objects\] batch=(\d+) site=(\w+) [^\n]*? hist=([0-9,]+) placed="
    for batch, site, values in re.findall(pattern, text):
        key = (batch, site)
        if key in seen or site not in result:
            raise ValueError("Duplicate or unknown profile row")
        seen.add(key)
        counts = list(map(int, values.split(",")))
        if len(counts) != 6:
            raise ValueError("Wrong histogram schema")
        result[site] = [a + b for a, b in zip(result[site], counts)]
    if not all(sum(result[SITES[i]]) > 0 for i in TUNABLE):
        raise ValueError("Missing training profile")
    return result


def quantile_cap(counts, quantile):
    threshold, seen = sum(counts) * quantile, 0
    for bound, count in zip(BOUNDS, counts):
        seen += count
        if seen >= threshold:
            return bound
    # The open-ended bin must not create an unbounded LDM allocation.
    return BOUNDS[-1]


def candidates(counts, seed, phase):
    result = {"seed": seed}
    for site in TUNABLE:
        if phase == "admit":
            values = {quantile_cap(counts[SITES[site]], .99)}
        elif phase == "caps":
            values = {0, quantile_cap(counts[SITES[site]], .90),
                      quantile_cap(counts[SITES[site]], .99)}
        else:
            values = {0, 4096, 8192}
        field = "reserves" if phase == "reserves" else "caps"
        for value in sorted(values):
            policy = seed.changed(field, site, value)
            if policy not in result.values():
                result[f"{SITES[site]}.{field}.{value}"] = policy
    return result


def make_plan(profile, seed=Policy(), phase="caps", repeats=2):
    if repeats < 2:
        raise ValueError("At least two independent screening blocks are required")
    counts = histogram(profile.read_text())
    policies = candidates(counts, seed, phase)
    blocks = []
    for repeat in range(repeats):
        names = [name for name in policies if name != "seed"]
        random.Random(20261001 + repeat).shuffle(names)
        # Three incumbent controls on each node expose within-block drift.
        names.insert(len(names) // 2, "seed")
        blocks.append(["seed"] + names + ["seed"])
    return {"schema": 1, "training_dataset": "train79", "phase": phase,
            "profile_sha256": hashlib.sha256(profile.read_bytes()).hexdigest(),
            "site_order": SITES, "payload_bytes": PAYLOAD,
            "policy_scope": "ordinary allocation admission; fixed audited scratch hints",
            "policies": {name: dataclasses.asdict(p) for name, p in policies.items()},
            "blocks": blocks}


def select(plan, rows, min_gain=.01):
    if not 0 < min_gain < 1:
        raise ValueError("Invalid minimum gain")
    rows = sorted(rows, key=lambda r: (r["block"], r["order"]))
    expected = {(b, order): name for b, names in enumerate(plan["blocks"])
                for order, name in enumerate(names)}
    indexed = {}
    for row in rows:
        key = (row["block"], row["order"])
        if key in indexed or expected.get(key) != row["policy"]:
            raise ValueError("Unexpected, duplicate or mislabelled measurement")
        if row["dataset"] != plan["training_dataset"]:
            raise ValueError("Validation data must not affect selection")
        if row.get("status") != "PASS" or row.get("profile") != 0:
            raise ValueError("Only successful unprofiled runs can select policies")
        if not math.isfinite(row["stage2"]) or row["stage2"] <= 0:
            raise ValueError("Invalid Stage2 measurement")
        digest = row.get("digest", [])
        if (len(digest) != 4 or not all(isinstance(x, str) for x in digest) or
                not all(re.fullmatch(r"[0-9]+", x) for x in digest[:2]) or
                not all(re.fullmatch(r"[0-9a-f]+", x) for x in digest[2:])):
            raise ValueError("Missing complete output fingerprint")
        indexed[key] = row
    if set(indexed) != set(expected):
        raise ValueError("Incomplete experiment; refusing partial winner selection")
    if len({tuple(r["digest"]) for r in rows}) != 1:
        raise ValueError("Correctness mismatch")
    # Do not pool absolute times across nodes; normalize to local controls.
    controls = {}
    noise = []
    for block in range(len(plan["blocks"])):
        nodes = {r.get("node") for r in rows if r["block"] == block}
        if len(nodes) != 1 or None in nodes:
            raise ValueError("A screening block must stay on one recorded node")
        times = [r["stage2"] for r in rows if r["block"] == block and r["policy"] == "seed"]
        if len(times) < 3:
            raise ValueError("Require start/middle/end controls in every block")
        median = statistics.median(times)
        controls[block] = median
        # Include full control range, not just MAD (three samples are sparse).
        noise.append((max(times) - min(times)) / median)
    threshold = max(min_gain, max(noise))
    scores = []
    for name in plan["policies"]:
        ratios = [r["stage2"] / controls[r["block"]] for r in rows if r["policy"] == name]
        gain = 1 - statistics.median(ratios)
        scores.append({"policy": name, "ratios": ratios, "gain": gain,
                       "eligible": name != "seed" and gain >= threshold and max(ratios) < 1})
    eligible = sorted((s for s in scores if s["eligible"]), key=lambda s: (-s["gain"], s["policy"]))
    winner = eligible[0]["policy"] if eligible else "seed"
    # This is a screening decision, not a claim of statistical significance.
    return {"selected": winner, "policy": plan["policies"][winner],
            "threshold": threshold, "local_controls": controls,
            "scores": scores, "requires_independent_confirmation": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("profile", type=Path)
    plan.add_argument("--seed", type=Path)
    plan.add_argument("--phase", choices=("caps", "reserves", "admit"), default="caps")
    plan.add_argument("--out", type=Path, required=True)
    pick = commands.add_parser("select")
    pick.add_argument("plan", type=Path)
    pick.add_argument("measurements", type=Path)
    pick.add_argument("--out", type=Path, required=True)
    env = commands.add_parser("env")
    env.add_argument("selection", type=Path)
    args = parser.parse_args()
    if args.command == "plan":
        seed = from_dict(json.loads(args.seed.read_text())["policy"]) if args.seed else Policy()
        if args.phase == "admit" and not args.seed:
            seed = dataclasses.replace(seed, caps=(0,) * len(SITES))
        result = make_plan(args.profile, seed, args.phase)
    elif args.command == "select":
        result = select(json.loads(args.plan.read_text()), json.loads(args.measurements.read_text()))
    else:
        policy = from_dict(json.loads(args.selection.read_text())["policy"])
        for name, value in policy.env().items():
            print(f"export {name}={value}")
        return
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
