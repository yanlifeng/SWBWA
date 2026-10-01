# Optional offline LDM admission selection

Production **B** is the built-in `CPE_LDM_MODE=3` policy. **C** refers to this
optional offline search using the same pool, not a separate allocator or an
online learner. It never changes the production default automatically.

Requirements: Python 3.7+ on the analysis machine, no extra packages. The tool
does not run alignment, submit jobs or access remote machines. It only writes
the explicitly named plan/selection files or prints bounded environment values.

## Workflow

1. On a training input, collect a complete mode-3 log with `SWBWA_LDM_PROFILE=1`
   and `-v 4`. Histograms count allocation events, not cache misses or accesses.
2. Generate a small candidate plan. Ordinary caps/reserves may change; audited
   ownership, scratch hints, pool capacity and allocator indices stay fixed.
3. Execute the plan's ordered blocks with profiling off. Each block stays on
   one recorded node and contains three seed controls to expose drift. Validate
   full output against an independent known-good fingerprint for each run.
4. Select only from a complete, correctness-checked measurement set. A candidate
   must improve every observation and exceed both the 1% gate and the largest
   relative control range. This is screening, not statistical significance.
5. Freeze the selection before independent inputs/nodes. Keep B if it does not
   improve the confirmation results. Never retune on the validation set.

```bash
python3 tools/ldm_policy.py plan training.log --phase caps --out plan.json
python3 tools/ldm_policy.py select plan.json measurements.json --out selection.json
python3 tools/ldm_policy.py env selection.json
```

`--phase admit` starts with ordinary admission off unless `--seed selection.json`
is given; existing manual scratch remains enabled. `--phase reserves` tunes only
reserves. `--seed` can chain explicit offline rounds; it is not loaded by SWBWA.

## Measurement schema

The JSON file contains an array with one row for every position in `plan.blocks`:

```json
[
  {
    "block": 0,
    "order": 0,
    "policy": "seed",
    "dataset": "train79",
    "node": "node-a",
    "status": "PASS",
    "profile": 0,
    "stage2": 40.0,
    "digest": ["5000000", "2107655733", "0123456789abcdef", "fedcba9876543210"]
  }
]
```

This row is a schema example, not a valid full experiment or reference checksum.
Use the plan's `training_dataset` label (currently `train79`) consistently.
`digest` holds the complete-output calls, bytes, sum and xor as strings. The
caller is responsible for reference validation before setting PASS; equality
between candidates alone cannot prove correctness. Missing/duplicate rows,
profiled runs, differing fingerprints or mixed nodes within a block are rejected.
Do not mix cropped-hash/CPE-copy-elided runs with the complete-output measurements.

## Preserved C candidate

The previously screened chain+reference candidate can be reproduced explicitly:

```bash
export SWBWA_LDM_CAPS=0,4096,1024,0,0,0,0,0,0,0,0
export SWBWA_LDM_RESERVES=0,4096,4096,0,4096,0,0,0,0,0,0
```

All other B settings are unchanged. This candidate was slightly slower than B
on the tested PE inputs and is not the default. Remove these overrides to return
to B. No online access profiling, live-object eviction or automatic escape
analysis is claimed; scripts and checks are reusable, the policy is bounded.
