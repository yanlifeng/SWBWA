# CPE LDM configuration

Use one build parameter: **`CPE_LDM_MODE=0|1|2`**, default **2**.

| Mode | Historical manual placement | New malloc pool | Intended use |
| --- | --- | --- | --- |
| 0 | Off; original heap fallbacks | Off | Explicit-LDM ablation |
| 1 | Off | 32 KiB per CPE, lifetime/access-tier policy | Low-intrusion allocation experiment |
| 2 | On | Off | Validated manual workspace optimizations |

This controls explicitly allocated private workspaces, not all physical LDM.
Runtime cache, the stack (`bsub -b`), static/compiler storage and registers do
not change. SIMD, scoring, scheduling and heap-workspace reuse also do not
change. The legacy worker2 raw-LDM path remains disabled in all three modes.

Mode 1 requires **non-MPI cgs_cross+pool**. Modes 0/2 preserve the existing
execution-mode restrictions of individual kernels. Mode 2 allows manual
placement; it does not override a separately disabled kernel feature. Keep
`CPE_KERNEL_OPT=1` for the complete optimized non-MPI cross+pool configuration.

```bash
bash build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 CPE_LDM_MODE=2
# Change only CPE_LDM_MODE to 0 or 1, rebuilding both passes each time.
make print-config EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=0 CPE_LDM_MODE=1
```

For matched Stage2 tests retain cache128, default chunk size and `-1`; use
`CPE_PROFILE=0 OUTPUT_MODE=discard CPE_DISCARD_DIGEST=0 DISCARD_HASH_BYTES=0`
and `SWBWA_DISCARD_HASH=1`. Full SAM construction/copy and fingerprints remain.

## Configuration ownership

- `include/swbwa_ldm_config.h`: mode, derived flags, capacity and budget checks.
- `src/slave/malloc_wrap.c`: single accounted SDK allocation gate. Mode 0
  refuses every site; mode 1 permits only the automatic arena; mode 2 retains
  manual calls. Every manual caller has a heap fallback.
- `src/slave/malloc_wrap.h`: wrapped malloc/calloc/realloc/free routing.
- `include/swbwa_ldm_alloc.h`: audited site and priority rules.
- `src/slave/ldm_alloc.c`: bitmap, ownership, realloc and arena lifetime.
- `src/slave/bwamem.c`, `slave.c`, `ksw.c`: existing lifecycle integration.

`SWBWA_CPE_MANUAL_LDM` and `SWBWA_CPE_LDM_ALLOC` are derived implementation
flags, not independent user options. Policy 4 retains the tested tiered
allocator unchanged. Historical profile/whitelist/full branches remain covered
by the isolated allocator harness, not exposed as production build modes.
Make rejects retired `CPE_LDM_ALLOC`, `CPE_MANUAL_LDM`, `CPE_LDM_BYTES` and
`LDM_SCRATCH_BUDGET` variables, including environment values, instead of
silently ignoring an old experiment command. Direct legacy -D flags also fail.

## Mode 1: bounded pool

The algorithm's allocation calls are unchanged. The wrapper uses `__func__`
to classify audited sites; constant comparisons can be folded by the compiler.
This is a placement policy, not automatic escape/access analysis.

| Priority | Objects |
| --- | --- |
| 1 | Temporary mate-rescue query profiles/DP; extension EH/QP; small global-DP scratch; context metadata; dedup sort/reorder scratch |
| 2 | SMEM intervals; chain seeds; chain sort/reversal scratch; reference slices and reg2aln query copies |
| Heap | SAM, CIGAR output, public/cached query profiles, unaudited sites |

A multiword bitmap manages 64-byte units. Requested lengths and spans support
ownership checks and bounded realloc copies. Free immediately reuses slots.
Realloc can move a growing LDM object to the original cross-segment segment-tree
pool; existing heap objects are never promoted. Above the old pool's 256 KiB
class limit, its unchanged system-allocation fallback remains in use.

Requests above half the payload go to heap. Priority-2 requests are limited
to 1/8 of the payload and must leave 1/4 free. Global-DP requests above 4 KiB
follow the lower-priority limits. No live pointer is evicted. Enough total
free bytes does not guarantee contiguous space. SDK refusal also falls back.

Payload, metadata and alignment slack share one **40 KiB** tracked request
budget. This is not total physical LDM use: stack, runtime storage and SDK
overhead require additional headroom.

### Lifetime

```text
alignment entry: begin arena, create context
  process multiple reads/pairs
  destroy context and reset reusable EH/QP pointers
  assert no live pool objects, suspend arena
SAM-copy entry on the same CPE: resume, end arena, release backing
```

SAM stays in heap and does not occupy LDM between entries. The existing task
list preserves CPE affinity. Main-core buffers and final shared output are never
replaced by private LDM. Verbosity 4 reports mode/capacity, placements, refusals,
spills and held bytes. Placement ratios are allocation statistics, not cache
hits or execution-time percentages.

## Mode 2: manual placement

Retains SMEM scratch/chain arena phase reuse, worker metadata, dedup scratch,
temporary KSW query profiles and bounded extension EH/QP placement. No new
automatic arena is created. The same budget and heap fallbacks apply.
Individual caps retain existing values. DP recurrence, tie order and scoring
do not change.

## Safety limits

Pool size and budget have internal, #ifndef-guarded overrides for regression
and research, not routine Makefile options. Payload above 96 KiB requires
manual placement disabled and separate cache/stack validation; 200 KiB payload
needs roughly 220 KiB backing. The experimental CPE guard refuses large arenas
unless at least 32 KiB remains below the entry stack marker. It does not prove
deepest-stack safety.

In the current cross runtime, completion parks without returning and the next
entry replaces PC without restoring SP. Entry/park frames can accumulate.
First-batch allocation success or standard athread persistence probes are not
a guarantee across many cross launches. Large pools are not defaults.
Changing that runtime protocol is separate from this configuration cleanup.

## Validation

```bash
python3 tests/test_ldm_modes.py
python3 tests/test_ldm_allocator.py
bash tests/run_cpe_kernel_checks.sh
```

Mode tests check SDK gating and rejected configurations. The allocator harness
keeps historical-policy coverage in a test-only configuration: randomized
ownership/realloc, zero sizes, fallback, large copy lengths and actual global-DP
score/CIGAR comparisons under sanitizers. Native tests do not validate Sunway
relocation or stack headroom. Target runs need both build passes, full output
fingerprints and per-batch arena/held-byte checks. Fingerprints are regression
checks, not order-sensitive proofs; full SAM MD5 is a separate validation.

### Matched PE measurements (2026-09-22)

One node, one process, 384 CPE, cgs_cross+pool, cache128, `-1`, default input
chunk size and full-SAM discard hashing with `CPE_DISCARD_DIGEST=0`.
Each configuration ran twice, with the second round in reverse order.
The table reports the faster Stage2 time in seconds; no input/output time
is subtracted from or added to Stage2.

| PE dataset | Mode 0 | Mode 1 | Mode 2 |
| --- | ---: | ---: | ---: |
| small_SRR7963242, 150 bp, 10,000,000 pairs | 172.341 | 161.458 | 158.193 |
| SRR2496709, 100 bp, 12,607,411 pairs | 60.636 | 58.020 | 56.259 |
| ERR1203383, 75 bp, 15,274,680 pairs | 63.767 | 62.145 | 59.543 |

All 18 full-SAM fingerprints (calls, bytes, sum and xor) matched the existing
baselines. No outstanding LDM allocations remained at batch completion.
The maximum within-configuration Stage2 spread was 0.14%; two repetitions
are not a confidence interval or a guarantee across other systems/datasets.
Mode 1 reduced Stage2 time by 2.54%-6.31% versus mode 0, but was 2.06%-4.37%
slower than mode 2. Manual placement therefore remains the default.
