# One-node to two-node scaling

The comparison uses the slowest rank's accumulated stage2 time as the stage2
makespan. The baseline is `../mpi_fine2_hash_1n6_retry_20260825` with one node
and six MPI ranks; this directory contains the two-node, 12-rank run with six
ranks per node.

All six `calls/bytes/hash-sum/hash-XOR` fingerprints match exactly between the
two runs. Discard mode created no SAM files.

| Dataset | Mode | 1 node / 6 ranks (s) | 2 nodes / 12 ranks (s) | Speedup | Efficiency | 6-rank max/min | 12-rank max/min |
|---|---:|---:|---:|---:|---:|---:|---:|
| ERR1203383 | PE | 68.403 | 37.013 | 1.848x | 92.4% | 1.056 | 1.147 |
| ERR1203383 | SE | 26.755 | 14.101 | 1.897x | 94.9% | 1.040 | 1.129 |
| SRR2496709 | PE | 70.547 | 37.374 | 1.888x | 94.4% | 1.115 | 1.222 |
| SRR2496709 | SE | 27.411 | 15.334 | 1.788x | 89.4% | 1.069 | 1.233 |
| small_SRR7963242 | PE | 232.978 | 118.240 | 1.970x | 98.5% | 1.024 | 1.039 |
| small_SRR7963242 | SE | 54.293 | 28.607 | 1.898x | 94.9% | 1.015 | 1.100 |

Mean stage2 speedup is **1.881x**, corresponding to **94.1% mean parallel
efficiency**. The larger 12-rank max/min ratios mainly reflect that useful work
per rank is roughly halved while the final indivisible expensive chunk remains
about the same size. Scheduler RMA remains small in absolute terms and is not
the dominant scaling limit in these runs.
