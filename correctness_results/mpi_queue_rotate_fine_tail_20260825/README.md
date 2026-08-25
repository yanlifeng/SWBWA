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
| ERR1203383 | PE | 1 | 42-48 | 29.450-32.004 | 63.061-66.903 | 1.061 | 60.339-64.106 | 0.997060-1.213642 | 255 / 4.993461 s | 269/271 | 38.890-47.533 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 2 | 40-47 | 28.039-29.751 | 62.867-66.471 | 1.057 | 60.070-63.913 | 0.916533-1.220743 | 255 / 5.013534 s | 266/271 | 43.710-49.363 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 3 | 40-48 | 29.066-31.922 | 62.686-66.440 | 1.060 | 59.994-63.831 | 0.906059-1.159854 | 255 / 5.013851 s | 268/271 | 43.604-46.699 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| SRR2496709 | PE | 1 | 41-50 | 27.321-31.174 | 61.363-65.414 | 1.066 | 58.763-62.565 | 0.998508-1.098337 | 100 / 9.412036 s | 273/274 | 42.325-52.250 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 2 | 42-48 | 25.848-30.957 | 61.257-67.358 | 1.100 | 58.392-64.880 | 0.943895-1.252972 | 100 / 9.459304 s | 270/274 | 44.613-48.224 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 3 | 40-52 | 22.609-27.707 | 60.373-67.889 | 1.124 | 57.331-65.373 | 0.997118-1.234463 | 100 / 9.479813 s | 268/274 | 47.580-55.193 | `56efd396c8f85a1e` | `9edd66de273acc48` |

## Main observations

- Repeated hash-enabled cases with matching aggregate fingerprints: 2/2.
- Largest scheduler-RMA/stage2 upper bound: 2.05% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 6/6. Run-level ranges and hashes remain complete.
- Slowest observed chunk: SRR2496709 PE chunk 100 at 9.479813 s.
- Lowest measured pipeline fread bandwidth: 38.890 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
