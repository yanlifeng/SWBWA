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
| ERR1203383 | PE | 1 | 40-47 | 27.735-30.221 | 64.643-68.147 | 1.054 | 62.008-65.703 | 0.847613-1.226167 | 255 / 4.995228 s | 269/271 | 45.131-50.311 | `0000000000000000` | `0000000000000000` |
| SRR2496709 | PE | 1 | 43-51 | 28.819-34.807 | 62.400-69.932 | 1.121 | 59.499-67.419 | 0.902282-1.141329 | 100 / 9.462080 s | 272/274 | 36.834-44.991 | `0000000000000000` | `0000000000000000` |

## Main observations

- Largest scheduler-RMA/stage2 upper bound: 1.90% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 2/2. Run-level ranges and hashes remain complete.
- Slowest observed chunk: SRR2496709 PE chunk 100 at 9.462080 s.
- Lowest measured pipeline fread bandwidth: 36.834 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
