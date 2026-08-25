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
| ERR1203383 | PE | 1 | 37-40 | 27.618-30.299 | 56.353-69.967 | 1.242 | 53.841-67.519 | 0.785014-0.981573 | 229 / 12.186040 s | 229/232 | 46.670-49.922 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 2 | 37-41 | 27.483-29.062 | 56.505-69.905 | 1.237 | 53.993-67.434 | 0.798257-1.158450 | 229 / 12.178719 s | 228/232 | 45.064-50.537 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 3 | 37-41 | 24.846-26.866 | 56.332-69.885 | 1.241 | 53.882-67.450 | 0.828406-0.994728 | 229 / 12.200566 s | 229/232 | 49.995-53.540 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| SRR2496709 | PE | 1 | 31-43 | 28.674-33.504 | 58.176-65.214 | 1.121 | 55.372-62.961 | 0.712632-1.106352 | 100 / 9.405718 s | 232/237 | 42.978-49.039 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 2 | 35-44 | 23.999-31.040 | 60.152-67.745 | 1.126 | 57.584-65.204 | 0.807427-0.973279 | 100 / 9.470594 s | 232/237 | 40.420-51.371 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 3 | 32-42 | 30.013-33.431 | 60.326-65.283 | 1.082 | 57.667-62.832 | 0.783073-1.043565 | 100 / 9.442312 s | 233/237 | 40.354-45.138 | `56efd396c8f85a1e` | `9edd66de273acc48` |

## Main observations

- Repeated hash-enabled cases with matching aggregate fingerprints: 2/2.
- Largest scheduler-RMA/stage2 upper bound: 2.05% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 6/6. Run-level ranges and hashes remain complete.
- Slowest observed chunk: ERR1203383 PE chunk 229 at 12.200566 s.
- Lowest measured pipeline fread bandwidth: 40.354 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
