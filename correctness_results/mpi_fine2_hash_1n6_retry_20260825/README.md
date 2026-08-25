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
| ERR1203383 | PE | 1 | 40-47 | 36.989-39.411 | 64.751-68.403 | 1.056 | 62.081-65.959 | 0.895981-1.335508 | 255 / 5.005719 s | 268/271 | 35.196-40.827 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | SE | 1 | 41-48 | 19.078-20.125 | 25.722-26.755 | 1.040 | 23.570-24.537 | 0.971336-1.259206 | 5 / 1.319505 s | 270/271 | 33.214-43.032 | `c4f605909a9eaa84` | `dd9be621d24957ae` |
| SRR2496709 | PE | 1 | 41-52 | 28.213-33.921 | 63.246-70.547 | 1.115 | 60.331-67.970 | 0.984358-1.181574 | 100 / 9.438087 s | 269/274 | 39.002-48.268 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | SE | 1 | 43-47 | 15.080-17.283 | 25.650-27.411 | 1.069 | 23.421-25.051 | 1.005531-1.114806 | 100 / 2.257153 s | 268/274 | 40.904-48.814 | `db44bb6c112040ba` | `337219a8b7ce62fa` |
| small_SRR7963242 | PE | 1 | 46-52 | 28.704-34.643 | 227.504-232.978 | 1.024 | 224.932-230.339 | 0.902225-1.318340 | 0 / 10.317269 s | 289/290 | 43.150-52.745 | `e31564c82ea77aa3` | `39ae379a71398f41` |
| small_SRR7963242 | SE | 1 | 47-50 | 11.583-13.361 | 53.483-54.293 | 1.015 | 51.127-51.889 | 1.116447-1.244836 | 0 / 3.042273 s | 286/290 | 57.928-67.072 | `9c841b9019c7d950` | `903fd6c660f2da8e` |

## Main observations

- Largest scheduler-RMA/stage2 upper bound: 4.90% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 6/6. Run-level ranges and hashes remain complete.
- Slowest observed chunk: small_SRR7963242 PE chunk 0 at 10.317269 s.
- Lowest measured pipeline fread bandwidth: 33.214 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
