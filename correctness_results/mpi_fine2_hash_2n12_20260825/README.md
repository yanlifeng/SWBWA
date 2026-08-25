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
| ERR1203383 | PE | 1 | 20-29 | 19.606-24.163 | 32.280-37.013 | 1.147 | 30.527-35.511 | 0.516246-0.988595 | 291 / 5.016644 s | 303/307 | 29.194-37.545 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | SE | 1 | 22-30 | 8.986-10.420 | 12.485-14.101 | 1.129 | 11.072-12.626 | 0.595362-1.109092 | 6 / 1.362477 s | 302/307 | 31.828-43.834 | `c4f605909a9eaa84` | `dd9be621d24957ae` |
| SRR2496709 | PE | 1 | 21-32 | 13.421-20.018 | 30.581-37.374 | 1.222 | 28.851-35.774 | 0.560795-1.247459 | 100 / 9.506673 s | 302/310 | 32.696-51.585 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | SE | 1 | 22-28 | 7.543-10.829 | 12.438-15.334 | 1.233 | 11.005-13.841 | 0.441462-0.862835 | 100 / 2.234643 s | 303/310 | 28.588-53.944 | `db44bb6c112040ba` | `337219a8b7ce62fa` |
| small_SRR7963242 | PE | 1 | 23-31 | 15.246-19.248 | 113.816-118.240 | 1.039 | 112.077-116.465 | 0.786240-1.312406 | 0 / 10.218810 s | 321/326 | 37.897-51.712 | `e31564c82ea77aa3` | `39ae379a71398f41` |
| small_SRR7963242 | SE | 1 | 23-31 | 6.207-8.781 | 26.013-28.607 | 1.100 | 24.652-26.996 | 0.702731-1.585782 | 0 / 2.967740 s | 322/326 | 43.500-82.313 | `9c841b9019c7d950` | `903fd6c660f2da8e` |

## Main observations

- Largest scheduler-RMA/stage2 upper bound: 8.88% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 6/6. Run-level ranges and hashes remain complete.
- Slowest observed chunk: small_SRR7963242 PE chunk 0 at 10.218810 s.
- Lowest measured pipeline fread bandwidth: 28.588 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
