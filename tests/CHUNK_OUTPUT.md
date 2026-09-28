# Chunk output and host workers

`process(step=2)` submits one input batch to `swbwa_output_write_chunk()`.
For MPI `single_unordered`, a nonempty chunk receives exactly one RMA extent,
even when it is larger or smaller than `SWBWA_OUTPUT_BUFFER_BYTES`. No batch
is combined with its successor. The writer sums the CPE-produced SAM lengths,
packs the records in read order, reserves the complete byte count, then uses `pwrite`
within that extent. Supplementary records stay with their source read.

The reusable packing allocation starts at `SWBWA_OUTPUT_BUFFER_BYTES` and
grows to the largest complete output chunk. This trades extra peak memory
for a single contiguous packing pass. It is not a fixed 64 MiB memory cap.
The input chunk limit and existing CPE SAM-buffer bounds still apply.
Legacy streaming callers of `swbwa_output_write()` retain buffered behavior.
Split/non-MPI output flushes at chunk boundaries but makes no RMA calls.
The length table is copied into batch-owned storage before the next CPE batch
can reuse its scratch arrays. It lives until Stage 3 completes, including in
the overlapped pipeline. Callers without a length table use a `strlen` fallback.

At `-v 4`, each extent reports input chunk ID/range and output offset/length.
In `single_unordered`, **global output is still unordered**: offsets are reserved
in completion order. The boundary-preserving interface also supports the
experimental `single_ordered` protocol described in [ORDERED_OUTPUT.md](ORDERED_OUTPUT.md).
Its monotonic-input requirement does not apply to the unordered writer.

## Six-MPE Work Sharing

Full-chip `cgs` and `cgs_cross` default to `HOST_MPE_THREADS=6`; `single` defaults
to `1`. Explicit `HOST_MPE_THREADS=1` retains the one-MPE control for either
full-chip mode. Six-MPE builds also default to six input readers; use
`SWBWA_INPUT_READERS=0` at runtime to retain serial `fread` without disabling
parallel SAM preparation/packing. One-MPE builds default to serial `fread`.
For cross execution use the normal complete two-pass `build.sh` workflow.
The job must allocate all six CGs to the process. Five persistent pthreads
call the SDK's `SET_TO_SWCG()` and verify placement using `rcid`. Failed
placement is an error, not silently reported as six-MPE execution.

The caller and helpers split contiguous read ranges for:

- Stage 2 SAM slice length sums and pointer/terminator assignment;
- Stage 3 SAM length prefix sums and complete chunk packing.

Helpers neither allocate/free application buffers nor enter MPI/athread.
The reader holds the shared helper team until its positioned reads finish.
Short SAM layout/packing calls use `swbwa_host_workers_run_ready()`: they take
the team if available, or execute all logical slices on the caller when busy.
The fallback preserves exactly the same slice boundaries, including when a
length pass uses helpers and its packing pass does not. It does not access
another dispatch's task storage. This avoids blocking computation behind I/O.
Task metadata remains live until every helper completes. CPE computation and
SAM generation are unchanged. Native tests exercise dispatch without pretending
to validate Sunway CG placement.

## Pipeline Ownership And Waiting

The overlapped pipeline uses bounded single-producer/single-consumer queues
with pthread mutexes and condition variables. Input capacity is one batch for
MPI dynamic scheduling, four otherwise; output capacity is four batches.
Closing a queue and publishing its final item use the same mutex. The consumer
drains all published items before observing EOF.

The processor waits for output queue space **before** generating a batch's
SAM. At that point at most three batches are queued and one is being written;
the new batch fits in the fifth SAM ring slot. Waiting only after generation
would permit overwriting a slot still owned by the writer. The writer releases
each batch after output completes, and performs the final flush once at EOF.

At `-v 4`, reader backpressure, processor input starvation, processor SAM-slot
waiting, and writer idle time are reported outside the stage callback timers.
These waits overlap across threads and must not be added together.
Direct pthread references are intentional: some Sunway static libstdc++ links
leave the weak `pthread_cond_broadcast` reference from `std::condition_variable`
unresolved. The queue must not depend on an unrelated translation unit pulling
that symbol into the executable.

## Tests

- `test_mpe_binding.c`: standalone Sunway six-CG placement, shared arrays,
  TLS/stack, synchronization, and join probe. Build on the login node; run only
  as an exclusively allocated one-node compute job.
- `test_host_workers.c`: concurrent dispatch from three callers (reader,
  Stage 2 preparation, writer), busy-team fallback, and repeated
  initialize/destroy cycles; included in `run_host_checks.sh`.
- `test_pipeline_queue.cpp`: empty input, FIFO order, slow reader/writer,
  one final flush, repeated pipelines, and SAM-ring ownership across wraparound;
  tested with both input capacities under ASan/UBSan.
- `python3 -B -m unittest discover -s tests -p test_output_rma.py -v`:
  three real MPI ranks, including an empty rank, real-file/RMA-only output,
  nonmonotonic chunk IDs, oversized chunks, within-chunk byte order,
  disjoint exhaustive extents, sample hashes, and one/six host threads.
- The same `test_output_rma.c` harness, built with `TEST_CHUNK_OUTPUT=1`,
  `SWBWA_HOST_MPE_THREADS=6` and real-file output, can run as a one-rank Sunway
  compute job. It independently checks all 2,883 output bytes, including every
  helper's slice. Its small synthetic payload is a packing test, not an
  end-to-end biological alignment accuracy test.

`OUTPUT_RMA_ONLY=1` still generates and packs complete SAM and performs RMA;
it omits disk writes and validates only the first 100 reads per chunk.
That sample check is not equivalent to full SAM accuracy verification.
