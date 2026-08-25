# MPI discard stage2 profile

Fixed configuration: `MPI_EXACT_READ_INDEX=1`, `MPI_INPUT_MODE=dynamic`, `OUTPUT_MODE=discard`, `-1`, `-v 4`. Hashes are combined offline across ranks with modulo-2^64 sum and XOR.

## Measurement semantics

- `stage2` and `CPE part3` are per-rank accumulated active time; they do not include time for another rank to finish.
- Stage1 bandwidth covers timed pipeline `fread` calls. The exact read-index scan and reference loading happen before the pipeline.
- Discard mode creates no SAM file and no output RMA window. With hashing enabled, each non-empty SAM blob passed to `swbwa_output_write` is hashed once.
- Per-rank hashes are printed independently. The checker combines them without MPI, using addition modulo 2^64 and XOR; both are independent of rank assignment and write order.
- With `-1`, stage3 hashing is outside stage2, although its cost still occurs before that rank asks the scheduler for its next chunk.
- Run-level ranges are collected independently of rank labels, so they remain valid when a batch system interleaves stderr from different nodes. `rank_summary.tsv` and `chunk_summary.tsv` are marked incomplete when that interleaving truncates or reorders detail rows.

## Run summary

| Dataset | Mode | Run | Chunks/rank | Stage1 (s) | Stage2 (s) | max/min | CPE part3 (s) | RMA (s) | Slowest chunk | Chunk rows | fread MiB/s | Hash sum | Hash XOR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| ERR1203383 | PE | 1 | 40-45 | 25.909-26.925 | 65.139-68.276 | 1.048 | 62.344-65.196 | 0.934275-1.070516 | 237 / 5.045165 s | 251/253 | 47.565-56.144 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 2 | 39-45 | 29.196-31.583 | 63.941-69.589 | 1.088 | 61.185-66.803 | 0.965755-1.176216 | 237 / 5.041897 s | 250/253 | 42.151-46.082 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 3 | 39-45 | 27.407-29.489 | 63.570-65.487 | 1.030 | 60.985-62.964 | 0.873701-1.185814 | 237 / 4.988705 s | 248/253 | 45.609-50.547 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| SRR2496709 | PE | 1 | 39-50 | 27.035-31.626 | 60.910-66.759 | 1.096 | 58.244-64.375 | 0.884335-1.067270 | 100 / 9.396935 s | 251/256 | 40.805-48.022 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 2 | 40-49 | 28.886-31.248 | 60.991-64.836 | 1.063 | 57.823-62.313 | 0.845225-1.100481 | 100 / 9.392007 s | 255/256 | 41.393-46.688 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 3 | 36-47 | 27.535-31.181 | 60.697-65.975 | 1.087 | 57.921-63.578 | 0.926701-1.115902 | 100 / 9.405283 s | 252/256 | 43.171-49.453 | `56efd396c8f85a1e` | `9edd66de273acc48` |

## Main observations

- Repeated hash-enabled cases with matching aggregate fingerprints: 2/2.
- Largest scheduler-RMA/stage2 upper bound: 1.87% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 6/6. Run-level ranges and hashes remain complete.
- Slowest observed chunk: SRR2496709 PE chunk 100 at 9.405283 s.
- Lowest measured pipeline fread bandwidth: 40.805 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
