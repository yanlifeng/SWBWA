# SWBWA MPI dynamic + single_unordered RMA performance study

Date: 2026-08-14

## 结论摘要

- 三种实现、has1/no1、PE/SE 共 12 项结果全部通过 MD5。
- 不做 `MPI_Iprobe` 时，PE 总时间为 123-132 秒；CPE 运行期间持续
  进入 MPI 后降到约 84 秒，说明 Sunway 上 RMA 的确依赖目标进程进入
  MPI 才能及时推进。
- `SERIALIZED + 用户锁` 与 `MULTIPLE + 无用户锁` 的四项总时间分别为
  227.114 秒和 227.945 秒，只差 0.37%，没有可辨认的性能胜负。
- `MULTIPLE` 将 PE no1 的输出 offset 申请从平均 3.947 秒降到 0.471
  秒，但总时间没有继续下降。等待转移到了 MPI 内部的 `MPI_Iprobe` 和
  `MPI_File_write_at`，最长单次 `MPI_Iprobe` 仍达到 9.576 秒。
- 如果保留当前三线程 pipeline，建议把 `probe + MPI_THREAD_MULTIPLE`
  作为工程基线，因为并发 MPI 调用语义正确且不需要全局用户锁，但不要
  把它视为性能优化。下一步应优先修正失效的 1 ms sleep，并测试较小的
  MPI-IO 写块或独立 MPI 服务线程。

## 1. Objective

This experiment compares three implementations of the CPE wait path while
keeping the alignment, dynamic input scheduler, and single-file unordered
output algorithms unchanged:

| Variant | CPE execution | MPI thread level | Process-local MPI mutex |
| --- | --- | --- | --- |
| `no_probe_serialized` | synchronous `swbwa_cpe_run()` | `MPI_THREAD_SERIALIZED` | yes |
| `probe_serialized` | asynchronous spawn, `MPI_Iprobe()` while CPE runs | `MPI_THREAD_SERIALIZED` | yes |
| `probe_multiple` | asynchronous spawn, `MPI_Iprobe()` while CPE runs | `MPI_THREAD_MULTIPLE` | no |

The `probe_*` variants execute the following logical loop:

```c
swbwa_cpe_spawn(...);
while (!cpe_finished) {
    MPI_Iprobe(...);
    usleep(1000);
}
swbwa_cpe_join();
```

The two serialized variants wrap MPI calls that may overlap across pipeline
threads, including `MPI_File_write_at`, with one process-local pthread mutex.
The MULTIPLE variant makes the same calls without that application mutex.

## 2. Test setup

- Platform: Sunway, one node, six MPI ranks, one CG per rank.
- Dataset: `SRR7963242`, tested as PE and SE.
- Input: MPI dynamic scheduler, exact read index enabled.
- Output: MPI `single_unordered` with a 64 MiB per-rank output buffer.
- Input chunk size: 19,660,800 FASTQ bytes, 49 chunks in total.
- Alignment thread count: one; all variants use the same CPE kernel.
- `has1`: command includes `-1`, therefore `kt_pipeline_single()` is used.
- `no1`: command omits `-1`, therefore the three-stage
  `kt_pipeline_queue()` is used.
- Timing reports are accumulated values. Pipeline stages overlap in `no1`,
  so stage times must not be added together.

Each source variant was built in an isolated directory. Jobs were submitted
sequentially, with at most one node occupied at any time.

## 3. Correctness

All 12 output checks passed after normalization and sorting:

| Variant | has1 PE | has1 SE | no1 PE | no1 SE |
| --- | ---: | ---: | ---: | ---: |
| `no_probe_serialized` | PASS | PASS | PASS | PASS |
| `probe_serialized` | PASS | PASS | PASS | PASS |
| `probe_multiple` | PASS | PASS | PASS | PASS |

Expected MD5 values:

- PE: `a799dc7268f389120ec92820e04b5118`
- SE: `0acfb1f46abd9fed3862c28a33e26da4`

The detailed results are in each variant's `has1/md5_results.tsv` and
`no1/md5_results.tsv` file. Validated SAM files were deleted; the complete
result directory is only about 1.1 MiB.

## 4. End-to-end results

Each cell contains `pipeline total / stage2 rank range`, in seconds.

| Variant | PE has1 | PE no1 | SE has1 | SE no1 |
| --- | ---: | ---: | ---: | ---: |
| `no_probe_serialized` | 132.298 / 43.035-77.869 | 123.559 / 78.262-92.762 | 37.518 / 11.744-15.501 | 32.906 / 12.143-16.041 |
| `probe_serialized` | 84.858 / 59.385-60.420 | 83.764 / 72.329-78.973 | 30.114 / 13.175-13.922 | 28.378 / 20.950-22.630 |
| `probe_multiple` | 84.542 / 60.025-60.989 | 83.322 / 73.192-77.573 | 30.707 / 11.867-14.886 | 29.374 / 14.298-22.541 |

The sum of the four pipeline totals is:

| Variant | Sum (s) | Change from no-probe |
| --- | ---: | ---: |
| `no_probe_serialized` | 326.281 | baseline |
| `probe_serialized` | 227.114 | 30.4% lower, 1.437x speedup |
| `probe_multiple` | 227.945 | 30.1% lower, 1.431x speedup |

Important observations:

- Adding MPI progress during CPE execution removes about 40-48 seconds from
  each PE run.
- Serialized and MULTIPLE are effectively tied. Their aggregate difference
  is 0.831 seconds, or 0.37%, which is smaller than normal run-to-run noise.
- MULTIPLE is 0.3-0.4 seconds faster for PE in this run. Serialized is
  0.6-1.0 seconds faster for SE.
- The no-probe PE runs have much wider stage2 ranges because delayed ticket
  acquisition also changes which ranks receive work.

## 5. Dynamic input balance

All variants process the same 49 chunks and 2,500,000 logical records.

| Case | Chunks per rank | Records per rank |
| --- | ---: | ---: |
| no-probe, PE has1 | 6-12 | 312,585-576,902 |
| probe serialized, PE has1 | 8-9 | 415,562-420,396 |
| probe MULTIPLE, PE has1 | 8-9 | 415,562-420,396 |
| no-probe, PE no1 | 7-10 | 364,380-473,368 |
| probe serialized, PE no1 | 8-9 | 415,562-420,398 |
| probe MULTIPLE, PE no1 | 8-9 | 415,562-420,398 |

The scheduler policy itself is not producing the large imbalance. Without
progress, ranks reach and complete the RMA ticket operation at very different
times. Once target ranks enter MPI during CPE work, the same scheduler gives
the expected 8-9 chunk distribution.

## 6. PE RMA and MPI-IO breakdown

The values below are the mean accumulated time per rank. `max call` is the
slowest single `MPI_Iprobe` seen anywhere in the run.

### PE has1

| Variant | Input ticket RMA (s) | Output reservation (s) | `MPI_File_write_at` (s) | Iprobe total (s) | Iprobe max call (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| no-probe serialized | 11.859 | 43.054* | 11.243* | 0.000 | 0.000 |
| probe serialized | 0.674 | 1.495 | 11.002 | 7.889 | 0.025 |
| probe MULTIPLE | 1.024 | 0.749 | 12.490 | 8.443 | 0.050 |

`*` One no-probe rank's output report was damaged by cross-rank stderr
interleaving, so these two means use five complete reports. The maximum remote
reservation was 58.699 seconds.

### PE no1

| Variant | Input ticket RMA (s) | Output reservation (s) | `MPI_File_write_at` (s) | Iprobe total (s) | Iprobe max call (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| no-probe serialized | 64.731 | 52.887 | 27.725 | 0.000 | 0.000 |
| probe serialized | 4.584 | 3.947 | 18.918 | 27.130 | 12.079 |
| probe MULTIPLE | 5.690 | 0.471 | 20.993 | 26.313 | 9.576 |

The no-probe serialized input ticket time consists mainly of:

- 51.949 seconds waiting for the process-local MPI mutex;
- 12.781 seconds inside `MPI_Win_unlock`.

The corresponding output reservation spends 48.989 seconds in
`MPI_Win_flush`. This is direct evidence that both input and output RMA are
waiting for MPI progress, not merely paying microsecond-scale atomic costs.

With MULTIPLE, the output reservation mean drops to 0.471 seconds, but the
whole run is not faster than serialized. The time reappears in blocking
`MPI_Iprobe` calls and in the MPI file-write path.

## 7. SE RMA and MPI-IO breakdown

For the no1 SE runs, the per-rank means are:

| Variant | Input ticket RMA (s) | Output reservation (s) | `MPI_File_write_at` (s) | Iprobe total (s) | Iprobe max call (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| no-probe serialized | 18.386 | 7.669 | 8.601 | 0.000 | 0.000 |
| probe serialized | 5.414 | 2.721 | 8.945 | 11.636 | 8.221 |
| probe MULTIPLE | 7.397 | 2.950 | 10.232 | 8.062 | 8.100 |

SE performs only about 10,000-14,000 probes per rank, yet an individual
`MPI_Iprobe` can still block for more than eight seconds. This rules out
"too many probes" as the only cause of the long tail.

## 8. What the Iprobe statistics mean

| Case | Probes/rank, mean | Time in probes/rank, mean | Slowest call | Calls >=1 s, all ranks |
| --- | ---: | ---: | ---: | ---: |
| serialized PE has1 | 846,614 | 7.889 s | 0.025 s | 0 |
| MULTIPLE PE has1 | 830,744 | 8.443 s | 0.050 s | 0 |
| serialized PE no1 | 601,346 | 27.130 s | 12.079 s | 27 |
| MULTIPLE PE no1 | 576,307 | 26.313 s | 9.576 s | 25 |
| serialized SE no1 | 9,791 | 11.636 s | 8.221 s | 25 |
| MULTIPLE SE no1 | 11,628 | 8.062 s | 8.100 s | 10 |

Two separate effects are visible:

1. The intended 1 ms polling interval is not reliably taking effect. A PE
   stage lasting about 60 seconds cannot execute roughly 830,000 iterations
   if every `usleep(1000)` actually sleeps for 1 ms. The current code ignores
   interrupted sleeps and does not measure their real duration.
2. In the pipelined no1 mode, some `MPI_Iprobe` calls themselves block for
   8-12 seconds. This remains true after removing the application mutex.
   Therefore Sunway MPI is still serializing or waiting for progress inside
   the library. `MPI_THREAD_MULTIPLE` permits concurrent calls; it does not
   guarantee asynchronous RMA progress or nonblocking implementation internals.

## 9. Causal interpretation

### 9.1 No probe + SERIALIZED

During a synchronous CPE kernel, the stage2 host thread does not enter MPI.
In no1 mode, stage1 and stage3 still run concurrently. A stage3 thread can
hold the process-local MPI mutex across a multi-second
`MPI_File_write_at`. The scheduler then waits tens of seconds for that mutex.

Even in has1, where there is no same-process pipeline contention, remote
`MPI_Win_unlock` and `MPI_Win_flush` wait until their target rank returns to
MPI. This matches the earlier standalone RMA test showing target-side progress
is required on this Sunway MPI implementation.

### 9.2 Probe + SERIALIZED

Calling `MPI_Iprobe` while CPEs run gives every rank opportunities to service
incoming RMA traffic. Input ticket and output reservation times fall sharply,
and the PE pipeline total falls from 123-132 seconds to about 84 seconds.

The remaining issue is lock convoying. In no1 mode the reader, CPE waiter, and
writer all contend for one MPI mutex. If one thread enters a long MPI call,
the progress thread cannot call `MPI_Iprobe`.

### 9.3 Probe + MULTIPLE

Removing the application mutex makes output offset reservation much more
predictable. In PE no1, its mean falls from 3.947 to 0.471 seconds. However,
end-to-end time does not improve. Individual `MPI_Iprobe` calls still block
for seconds, and MPI file writes still consume 15-25 seconds per rank.

The user mutex was one bottleneck, but not the root implementation limit.
The Sunway MPI library's internal progress and serialization behavior remains
on the critical path.

## 10. Which version is best now?

There is no meaningful performance winner between the two probe variants:

- Purely by the sum of these four runs, `probe_serialized` is 0.37% faster.
- For both PE tests, `probe_multiple` is about 0.5% faster.
- The difference is too small for a one-run benchmark to distinguish from
  system noise.

For the current three-thread pipeline, `probe_multiple` is the better
engineering baseline because concurrent MPI calls are standards-compliant
without a global application mutex, and its output reservation path is much
less variable. It should not be described as a performance optimization on
this platform.

If a single output file is not mandatory, dynamic input plus split output is
still the practical high-throughput configuration because it removes this
shared RMA/MPI-IO progress problem entirely.

## 11. Recommended next steps

### Priority 1: make polling deterministic

Replace unchecked `usleep(1000)` with an EINTR-safe monotonic sleep, preferably
`clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, ...)`, and record actual sleep
time. This should eliminate hundreds of thousands of unnecessary MPI calls.
It will reduce pressure, but cannot by itself fix one `MPI_Iprobe` blocking for
8-12 seconds.

### Priority 2: bound blocking MPI-IO calls

The 64 MiB buffer produces a small number of very long
`MPI_File_write_at` calls. Test 4, 8, 16, and 64 MiB write units. Smaller units
increase call count but let the process re-enter progress code between writes.
The useful metric is total pipeline time plus the maximum single Iprobe/RMA
completion time, not only aggregate file bandwidth.

### Priority 3: create one MPI service thread

Make one host thread own input ticket RMA, output offset RMA, and progress.
Pipeline threads submit requests through local queues. This removes both the
application lock convoy and uncontrolled concurrent calls into Sunway MPI.

The service thread must not spend many seconds in one blocking file write.
Two candidates should be tested separately:

- reserve offsets in the MPI service thread, then use POSIX `pwrite` from a
  writer thread to non-overlapping ranges, if the parallel filesystem supports
  this efficiently;
- issue smaller MPI file writes and interleave them with progress calls.

`MPI_File_iwrite_at` is worth a microbenchmark, but the existing RMA results
suggest it may not progress asynchronously without further MPI calls.

### Priority 4: isolate input and output in microbenchmarks

Run the same instrumented binary with:

- dynamic input + split output, to measure scheduler RMA alone;
- static input + single output, to measure offset RMA and MPI-IO alone;
- RMA-reserved offsets + POSIX `pwrite`, to separate MPI-IO from filesystem
  behavior.

This will show whether the next optimization should target MPI progress,
MPI-IO, or the filesystem write size.

## 12. Raw data

The raw logs are stored under:

```text
correctness_results/mpi_rma_singleout_study_20260814/
  no_probe_serialized/{has1,no1}/
  probe_serialized/{has1,no1}/
  probe_multiple/{has1,no1}/
```

Some stderr lines from different ranks are interleaved, but all six total and
stage2 timing rows and the dedicated RMA/progress counters used above are
present. Missing or damaged stage3 rows were not used in the conclusions.
