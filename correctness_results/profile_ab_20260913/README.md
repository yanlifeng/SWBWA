# SWBWA stage-2 kernel profiling: old vs new

Purpose: measure what the CPE optimizations between `fc658b8` and `a3e7beb` actually
bought in **stage 2**, and identify where stage-2 time goes **now**, so the next
optimization round can target real hotspots.

## What is compared

| | commit | date | role |
| --- | --- | --- | --- |
| `old` | `fc658b8` | 2026-08-24 | produced the historical baseline numbers |
| `new` | `a3e7beb` | 2026-09-09 | current `master` (`HEAD`) |

Both trees were built and run **on the same cluster, back-to-back, in the same
`bsub` window**, so node / clock / filesystem state cancels. Every condition is a
single blocking `bsub -I -b`; at most one node was ever in use.

Built for each variant (`JOBS=8`):
`{old,new} × {single(system), cgs_cross(pool)} × {CPE_PROFILE=0, CPE_PROFILE=1}`.
The `CPE_PROFILE=1` builds only add LWPF counters; they do not change alignment
results (verified by MD5, see below).

Each run writes SAM to the shared filesystem and deletes it immediately. Only logs
were transferred back (4.7 MB total, throttled to 500 KB/s).

- `has1` (`-1`) everywhere, so stage 1/3 hold all the IO and stage 2 is pure CPE compute.
- Datasets: `ERR1203383`, `SRR2496709`, `small_SRR7963242`; both SE and PE.
- Condition order per condition: `old_np → new_np → new_p → old_p`, where
  `np` = clean stage-2 wall clock, `p` = LWPF kernel counters. The ABBA ordering
  cancels monotonic drift.

## Result 1 — stage 2 wall clock (the headline)

| config | mean delta | range |
| --- | ---: | --- |
| `cgs_cross + pool` | **−12.8%** | −7.1% … −21.9% |
| `single + system` | **−12.7%** | −5.7% … −24.1% |

Per-dataset detail is in `report_cross.md` / `report_single.md`. The two configs
agree closely, which is a good sign that the improvements are not config-specific.

## Result 2 — kernel cycles, old vs new

`old` instruments **9 regions**; `new` instruments **33** (the original 9 plus 24
nested sub-regions). The first 9 are a prefix of the new list with identical names
and order (`fc658b8:slave/swbwa_cpe_profile.h` vs
`HEAD:src/slave/swbwa_cpe_profile.h`), so those 9 are directly comparable.

Regions that clearly improved:

| region | cross | single |
| --- | --- | --- |
| `MATE_RESCUE` (PE only) | −20% … −32% | −44% … −56% |
| `SAM_FORMAT` (PE) | −18% … −31% | −36% … −55% |
| `MEM_CHAIN` | −9% … −13% | −4% … −8% |
| `CHAIN_EXTENSION` | −5% … −23% | −6% … −9% |
| `PAIRING` (PE only) | −23% … −30% | −17% … −30% |

`SAM_COPY` looks like a large regression in cross (+53% … +127%), but it is
**re-scoping**: code moved between `SAM_FORMAT` and `SAM_COPY`. Merged they are:

| config | PE | SE |
| --- | ---: | ---: |
| cross | −17% … −31% | −2% … −5% |
| single | −36% … −55% | ~0% |

### Two cross-only regressions worth a look

| region | cross | single |
| --- | --- | --- |
| `ALIGNMENT_FINALIZE` | **+24% … +33%** (6/6 conditions) | −55% … −64% |
| `CHAIN_FILTER` | −8% … **+4%** | −11% … −24% |

`ALIGNMENT_FINALIZE` gets consistently *slower* in `cgs_cross` but dramatically
faster in `single`. Since the two builds differ in `EXEC_MODE` and allocator, this
looks like a cross-mode-specific effect rather than a general regression, and is
the most concrete lead for the next round.

## Result 3 — where the time goes now (new build, share of CPE total)

The LWPF regions are **hierarchical**, not a flat list: a region's counter includes
any region opened inside it. The tree below was read off the sources
(`src/slave/slave.c`, `src/slave/bwamem.c`, `src/slave/bwamem_pair.c`):

```text
WORKER_ALIGNMENT                    worker12_s_pre_fast  (whole pre-SAM phase)
├── MEM_CHAIN                       mem_align1_core_impl
│   ├── MEM_CHAIN_COLLECT           mem_chain -> mem_collect_intv
│   │   └── MEM_COLLECT_{FIRST,SPLIT,LAST,SORT}
│   └── MEM_CHAIN_BUILD             mem_chain
│       └── CHAIN_BUILD_{SA,RID,TREE_SEARCH,MERGE,INSERT,FINALIZE,REPETITIVE}
├── CHAIN_FILTER
├── CHAIN_EXTENSION
│   └── CHAIN_EXTENSION_DP          mem_chain2aln
├── ALIGNMENT_FINALIZE
└── SAM_FORMAT                      worker12_pre_fast -> mem_sam_pe
    ├── MATE_RESCUE                 mem_sam_pe (PE only)
    │   ├── MATE_REF_FETCH
    │   ├── MATE_KSW_ALIGN          swbwa_matesw_run_one/_pair
    │   │   └── KSW_QUERY_INIT_{FWD,REV}, KSW_DP_{FWD,REV}
    │   └── MATE_DEDUP
    │       └── DEDUP_{SORT_END,REDUNDANCY,SORT_SCORE}
    └── PAIRING                     mem_sam_pe (PE only)
SAM_COPY                            worker12_s_fast  (true sibling)
```

**`SAM_FORMAT` is a child of `WORKER_ALIGNMENT`** (opened and closed inside
`worker12_pre_fast`, which itself runs inside the `WORKER_ALIGNMENT` bracket), and
in PE mode it **contains** `MATE_RESCUE` and `PAIRING`. Only `SAM_COPY` is
disjoint from `WORKER_ALIGNMENT`, so the CPE total used here is
**`WORKER_ALIGNMENT + SAM_COPY`**. An earlier revision of this report wrongly used
`WORKER_ALIGNMENT + SAM_FORMAT + SAM_COPY` as the denominator, which double
counted `SAM_FORMAT`; the numbers below are corrected.
Verified numerically: `sum(WORKER_ALIGNMENT's children) / WORKER_ALIGNMENT` is
0.98–1.00 across all 12 conditions.

**PE, errors ~0 (ERR1203383) vs high-error (small_SRR7963242), % of CPE total:**

| region | ERR1203383 | SRR2496709 | small_SRR7963242 |
| --- | ---: | ---: | ---: |
| `WORKER_ALIGNMENT` | 98.8% | 98.9% | 99.7% |
| └ `MEM_CHAIN` | **46.6%** | **50.3%** | 25.9% |
| └ `CHAIN_EXTENSION` | 11.1% | 14.6% | 22.0% |
| └ `SAM_FORMAT` | **34.0%** | 29.2% | **45.7%** |
| &nbsp;&nbsp;└ `MATE_RESCUE` | 25.3% | 19.5% | 42.4% |
| &nbsp;&nbsp;&nbsp;&nbsp;└ `MATE_KSW_ALIGN` | 15.8% | 14.6% | 31.4% |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└ `KSW_DP_FORWARD` | 12.6% | 12.2% | 22.3% |
| &nbsp;&nbsp;└ `PAIRING` | 1.1% | 0.9% | 0.3% |
| `SAM_COPY` | 1.2% | 1.1% | 0.3% |

**SE, % of CPE total:**

| region | ERR1203383 | SRR2496709 | small_SRR7963242 |
| --- | ---: | ---: | ---: |
| `MEM_CHAIN` | **63.9%** | 63.4% | 49.6% |
| └ `MEM_CHAIN_COLLECT` | 37.6% | 38.4% | 23.8% |
| └ `MEM_CHAIN_BUILD` | 26.3% | 25.0% | 25.8% |
| `CHAIN_EXTENSION` | 14.2% | 16.5% | 34.1% |
| └ `CHAIN_EXTENSION_DP` | 5.7% | 9.6% | 25.8% |
| `SAM_FORMAT` | 11.0% | 10.7% | 5.4% |
| `MATE_RESCUE` | 0.0% | 0.0% | 0.0% |

### What `SAM_FORMAT` actually is

It is **not** just string formatting. In `worker12_pre_fast` (src/slave/bwamem.c)
the region brackets the call to `mem_sam_pe`, i.e. the whole per-pair paired-end
resolution:

1. mate rescue — SW re-alignment of the mate against candidate regions
   (`swbwa_mem_matesw_dual` / `_one`), which is the KSW DP work counted under
   `MATE_RESCUE` → `MATE_KSW_ALIGN` → `KSW_DP_*`;
2. dedup/sort of the rescued candidates (`MATE_DEDUP`);
3. pairing the two ends using the insert-size stats (`PAIRING`);
4. primary marking, mapQ computation, and finally emitting the SAM text for both
   mates.

So of PE's `SAM_FORMAT`, ~74% (ERR1203383) to ~93% (small_SRR7963242) is
`MATE_RESCUE`, i.e. **real alignment compute**, and only the remainder is
bookkeeping plus string building. In SE mode `SAM_FORMAT` merely wraps
`mem_mark_primary_se` + `mem_reg2sam`, which is why it is 5–11% there.

**Consequence:** `MEM_CHAIN` — not `SAM_FORMAT` — is the largest single region in
PE as well (26–50%). `SAM_FORMAT` is second, but most of it *is* the mate-rescue
KSW, which the report counts separately as `MATE_RESCUE`. Do not add
`SAM_FORMAT` + `MATE_RESCUE`: the latter is inside the former.

`small_SRR7963242` stands out again: its PE `SAM_FORMAT` is 45.7% (dominated by
mate rescue at 42.4%), and its SE `CHAIN_EXTENSION_DP` reaches 25.8% — it is also
the dataset that improved least, i.e. the current optimizations do not target
inner-extension DP, which this dataset stresses hardest.

## Recommended next targets

Ranges below span all 12 conditions, `% of CPE total = WORKER_ALIGNMENT + SAM_COPY`.

1. **`MEM_CHAIN`** — 26–67%, the largest region in both modes. Split into
   `MEM_CHAIN_COLLECT` (12–39%, SMEM backward search on the BWT) and
   `MEM_CHAIN_BUILD` (14–29%, of which `CHAIN_BUILD_SA` alone is 10–19%).
2. **`SAM_FORMAT` / `MATE_RESCUE` (PE only)** — 29–46% inclusive; the actionable
   part is `MATE_KSW_ALIGN` (15–32%) → `KSW_DP_FORWARD` (12–23%). Remember
   `MATE_RESCUE` is *inside* `SAM_FORMAT`.
3. **`CHAIN_EXTENSION_DP`** — 5–26%, peaking on `small_SRR7963242`; barely
   optimized so far.
4. **`CHAIN_BUILD_SA`** — 10–19%.
5. Investigate the **cross-only `ALIGNMENT_FINALIZE` / `CHAIN_FILTER` regressions.**

## Correctness

MD5 is over the header-free SAM, compared against
`scripts/bigdata_expected_md5.tsv`. Both the SWBWA run and the `md5sum` are
submitted as their own blocking `bsub` jobs (SWBWA must be the direct child
command), so nothing large runs on the login node. SAMs are deleted after checking.

Verified that the `CPE_PROFILE=1` builds do not change the result, using
`cross + pool`, `ERR1203383`, `PE` (expected `dc5c0a6babd41641db22808caedb7a44`):

| variant | SAM bytes | md5 | result |
| --- | ---: | --- | --- |
| `old_p` (`fc658b8`) | 8029871401 | `dc5c0a6babd41641db22808caedb7a44` | **PASS** |
| `new_p` (`a3e7beb`) | 8029871401 | `dc5c0a6babd41641db22808caedb7a44` | **PASS** |

Both variants also produced byte-identical SAM sizes in **every** condition of the
A/B matrix, which is a second independent consistency signal. Details in
`md5check.log`. The SAMs were removed afterwards; residual disk use on the cluster
is `runs/` 326 MB (binaries) + `logs/` 5.1 MB.

## Method note: why runs must use `bsub ... ./SWBWA` directly

On Sunway the CPE/mpe launcher only sets up the CPE mapping for the **direct child
process**. Launching SWBWA as `bsub ... bash script.sh` makes it a grandchild, and
it dies with SIGSEGV within a second (the baked-in layout addresses do not match).
This is what broke an earlier attempt at this experiment; it was not a code bug.

Also note: a bare `./SWBWA` with no arguments segfaults even in known-good builds,
so "it crashes with no args" is not evidence of a broken build.

## Files

| file | content |
| --- | --- |
| `report_cross.md` / `report_single.md` | rendered tables per config |
| `ab_stage2.tsv` | stage 1/2/3 + total, old vs new, all 12 conditions |
| `kernel_diff.tsv` | the 9 comparable regions, old/new cycles + delta + share |
| `kernel_new_breakdown.tsv` | all 33 regions of the new build, share of CPE total |
| `kernel_cycles.tsv` | raw per-region cycles, both variants, all conditions |
