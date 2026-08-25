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
| ERR1203383 | PE | 1 | 37-40 | 25.640-27.604 | 60.030-69.771 | 1.162 | 57.683-67.339 | 0.001294-0.210893 | 229 / 12.178161 s | 228/232 | 44.882-51.550 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 2 | 37-40 | 25.088-27.349 | 60.017-69.778 | 1.163 | 57.681-67.341 | 0.001303-0.199967 | 229 / 12.172337 s | 230/232 | 47.788-53.278 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | PE | 3 | 37-40 | 23.747-25.781 | 60.025-70.212 | 1.170 | 57.677-67.612 | 0.027274-0.135512 | 229 / 12.168871 s | 230/232 | 48.315-56.065 | `dae1072cac7f170d` | `1fe63f8e1e740aaf` |
| ERR1203383 | SE | 1 | 38-39 | 13.638-15.651 | 24.712-26.816 | 1.085 | 22.713-24.809 | 0.001803-0.183145 | 229 / 2.127569 s | 228/232 | 42.501-50.852 | `c4f605909a9eaa84` | `dd9be621d24957ae` |
| ERR1203383 | SE | 2 | 38-39 | 13.727-15.512 | 24.611-26.920 | 1.094 | 22.595-24.884 | 0.067483-0.141697 | 229 / 2.147241 s | 230/232 | 44.704-50.487 | `c4f605909a9eaa84` | `dd9be621d24957ae` |
| ERR1203383 | SE | 3 | 38-39 | 14.638-15.585 | 24.769-26.915 | 1.087 | 22.703-24.887 | 0.001687-0.226986 | 229 / 2.155756 s | 230/232 | 42.474-47.912 | `c4f605909a9eaa84` | `dd9be621d24957ae` |
| SRR2496709 | PE | 1 | 28-45 | 28.968-37.548 | 61.418-70.364 | 1.146 | 58.899-68.330 | 0.011289-0.319936 | 100 / 9.406111 s | 235/237 | 39.232-46.385 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 2 | 27-46 | 30.056-37.987 | 61.261-70.285 | 1.147 | 58.694-68.281 | 0.058301-0.359054 | 100 / 9.413651 s | 234/237 | 35.818-47.178 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | PE | 3 | 30-44 | 29.408-37.627 | 61.210-71.112 | 1.162 | 58.642-69.056 | 0.048699-0.293686 | 100 / 9.463655 s | 234/237 | 38.021-45.966 | `56efd396c8f85a1e` | `9edd66de273acc48` |
| SRR2496709 | SE | 1 | 39-40 | 13.104-16.801 | 24.625-28.253 | 1.147 | 22.642-26.191 | 0.002151-0.120571 | 100 / 2.228845 s | 232/237 | 40.822-49.760 | `db44bb6c112040ba` | `337219a8b7ce62fa` |
| SRR2496709 | SE | 2 | 38-41 | 13.976-17.982 | 24.378-28.217 | 1.157 | 22.391-26.101 | 0.023268-0.245264 | 100 / 2.230218 s | 233/237 | 37.910-47.803 | `db44bb6c112040ba` | `337219a8b7ce62fa` |
| SRR2496709 | SE | 3 | 38-41 | 14.716-18.232 | 24.673-28.188 | 1.142 | 22.637-26.157 | 0.069953-0.211811 | 100 / 2.231607 s | 232/237 | 36.463-44.304 | `db44bb6c112040ba` | `337219a8b7ce62fa` |
| small_SRR7963242 | PE | 1 | 42-43 | 18.411-25.145 | 227.157-235.433 | 1.036 | 224.632-232.773 | 0.001225-0.135136 | 0 / 10.354868 s | 250/254 | 55.032-69.288 | `e31564c82ea77aa3` | `39ae379a71398f41` |
| small_SRR7963242 | PE | 2 | 42-43 | 19.999-25.761 | 226.095-232.131 | 1.027 | 223.657-229.771 | 0.001190-0.164949 | 0 / 10.273796 s | 251/254 | 53.529-67.834 | `e31564c82ea77aa3` | `39ae379a71398f41` |
| small_SRR7963242 | PE | 3 | 42-43 | 21.877-25.988 | 226.709-231.034 | 1.019 | 224.242-228.470 | 0.001260-0.148024 | 0 / 10.249471 s | 250/254 | 53.200-67.082 | `e31564c82ea77aa3` | `39ae379a71398f41` |
| small_SRR7963242 | SE | 1 | 42-43 | 9.300-10.915 | 52.522-53.945 | 1.027 | 50.307-51.776 | 0.001200-0.191454 | 0 / 2.944834 s | 250/254 | 61.310-76.739 | `9c841b9019c7d950` | `903fd6c660f2da8e` |
| small_SRR7963242 | SE | 2 | 42-43 | 11.198-13.009 | 52.336-53.800 | 1.028 | 50.123-51.650 | 0.001937-0.192296 | 0 / 2.936674 s | 250/254 | 50.418-66.425 | `9c841b9019c7d950` | `903fd6c660f2da8e` |
| small_SRR7963242 | SE | 3 | 42-43 | 9.702-11.502 | 52.291-53.755 | 1.028 | 50.096-51.507 | 0.001161-0.190115 | 0 / 2.969355 s | 252/254 | 56.048-77.188 | `9c841b9019c7d950` | `903fd6c660f2da8e` |

## Main observations

- Repeated hash-enabled cases with matching aggregate fingerprints: 6/6.
- Largest scheduler-RMA/stage2 upper bound: 1.01% (max RMA divided by min stage2 within a run).
- Runs with incomplete rank/chunk detail parsing due to combined-log interleaving: 18/18. Run-level ranges and hashes remain complete.
- Slowest observed chunk: ERR1203383 PE chunk 229 at 12.178161 s.
- Lowest measured pipeline fread bandwidth: 35.818 MiB/s.

Detailed data are in `rank_summary.tsv` and `chunk_summary.tsv`; the original per-rank reports remain in the `.runNN.log` files.
