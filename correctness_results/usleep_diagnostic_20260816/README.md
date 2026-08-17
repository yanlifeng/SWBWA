# Sunway CPE progress sleep diagnostic

## Result

`usleep(1000)` itself is valid, but it is not reliable as a polling throttle
while the real alignment CPE kernel is active. Runtime signals interrupt most
calls and the existing code ignores `EINTR`, so it immediately executes the
next `MPI_Iprobe`.

In the 300,000-read SE diagnostic, the six ranks made 12,280 sleep calls during
alignment. Of those, 11,849 (96.5%) returned `-1` with `errno == EINTR`. The
shortest interrupted call returned after about 12 us. The same executable's
short format and SAM-generation kernels, where these signals were absent,
normally slept for roughly 6-10 ms.

An independent synthetic CPE test did not reproduce the signals: `usleep`,
relative `nanosleep`, and absolute `clock_nanosleep` all completed normally and
slept for about 10 ms. This also shows that the node's effective timer/scheduler
granularity is much coarser than the requested 1 ms.

The tests establish interruption as the cause, but do not identify the signal
number. Replacing runtime-installed signal handlers merely to observe it could
change athread behavior and was deliberately avoided.

## Comparison

All cases used one node, six MPI ranks, dynamic input, single unordered output,
one CG per rank, and the same 300,000-read `SRR7963242` SE subset.

| Wait method | Mean alignment probes/rank | Effective interval | Stage 2 range |
| --- | ---: | ---: | ---: |
| Original `usleep(1000)`, ignore `EINTR` | 2,047 | about 0.85-0.97 ms in this subset | 1.843-2.935 s |
| Relative `nanosleep`, retry kernel `remaining` | 942 | about 1.94-2.11 ms | 1.841-2.939 s |
| Absolute `clock_nanosleep` | 190 | about 9.99-10.21 ms | 1.844-2.931 s |
| Monotonic 1 ms spin with `sched_yield` | 1,902 | about 1.00 ms | 1.834-2.920 s |
| Monotonic deadline plus relative `nanosleep` | 1,006 | about 1.82-1.97 ms | 1.840-2.933 s |
| Previous relative-sleep helper | 1,014 (all CPE phases) | about 1.82-1.97 ms | 1.844-2.934 s |

The yield loop controls the interval accurately but performs roughly 1.3-2.4
million `sched_yield` calls per alignment kernel. Absolute sleeping is cheap
and polls approximately once per 10 ms on this platform. A full SWBWA run
below showed that this lower polling rate did not hurt dynamic-input or
single-file-output performance.

The selected implementation now uses an absolute `CLOCK_MONOTONIC` deadline
through `clock_nanosleep(..., TIMER_ABSTIME, ...)`. Interrupted calls retry the
same 1 ms deadline, so signals cannot either trigger a probe storm or extend
the intended wait. The effective wait remains about 10 ms because of the
platform's timer granularity.

## Full SWBWA verification

The absolute-sleep production build was also tested on the complete
`SRR7963242` dataset with six MPI ranks, dynamic input, exact read indexing,
and `single_unordered` output. Both PE and SE passed their expected MD5 in
has1 and no1 modes.

Each performance cell is `pipeline total / stage 2 rank range`, in seconds.
The old baseline is the `probe_multiple` build that used `usleep(1000)` and
ignored `EINTR`.

| Wait method | PE has1 | PE no1 | SE has1 | SE no1 |
| --- | ---: | ---: | ---: | ---: |
| Old `usleep` baseline | 84.542 / 60.025-60.989 | 83.322 / 73.192-77.573 | 30.707 / 11.867-14.886 | 29.374 / 14.298-22.541 |
| Absolute `clock_nanosleep` | 78.456 / 56.777-58.569 | 74.236 / 64.190-70.882 | 27.675 / 12.320-14.811 | 25.362 / 14.739-20.840 |

The absolute version reduced mean probes per rank from 830,744 to 3,369 in
PE has1 and from 576,307 to 2,356 in PE no1. It did not regress any pipeline
total or stage 2 range in this run. These are single-run results on different
nodes, so the apparent 7-14% total-time improvement should be treated as
encouraging rather than a precise speedup claim.

Long MPI calls remain in pipelined no1 mode. The slowest `MPI_Iprobe` was
9.240 seconds for PE and 7.383 seconds for SE even though probe counts were
low. This confirms that the remaining tail comes from Sunway MPI progress and
internal serialization, not from excessive polling or the sleep API.

## Logs

- `test_usleep_progress.log`: independent synthetic CPE test.
- `run.log`: instrumented original `usleep` behavior.
- `run_nanosleep2.log`: relative `nanosleep` using returned remaining time.
- `run_clock_nanosleep.log`: absolute `clock_nanosleep`.
- `run_yield_wait.log`: monotonic deadline with `sched_yield`.
- `run_deadline_nanosleep.log`: monotonic deadline with recomputed relative sleep.
- `run_fix_verify.log`: previous relative-sleep production helper.
- `../clock_nanosleep_full_20260816/`: full-dataset absolute-sleep logs and
  MD5 results.

Remote jobs: `7794790`, `7794818`, `7794840`, `7794855`, `7794872`, and
`7794897`. Full-dataset run jobs: `7795298`, `7795322`, `7795333`, and
`7795350`.
