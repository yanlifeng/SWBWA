# Experimental Ordered MPI Output

Build with `USE_MPI=1 MPI_INPUT_MODE=dynamic OUTPUT_MODE=single_ordered`.
For cross execution use the complete two-pass `build.sh` workflow. The default
output mode is unchanged. Static input is rejected for this first implementation.
"Ordered" means FASTQ input/read order, not genomic-coordinate sorting.

## Protocol

- Global input tickets assign monotonically increasing chunk IDs to each rank.
- Rank 0 exposes `chunk_count + 1` 64-bit prefix slots. Slot 0 starts at zero;
  the others start at `UINT64_MAX` (not ready).
- A writer measures/packs its chunk, waits for its prefix with atomic
  `MPI_Fetch_and_op(MPI_NO_OP)` plus `MPI_Win_flush`, and publishes the next
  offset with `MPI_Accumulate(MPI_REPLACE)` plus flush.
- Publishing the next offset does not wait for disk I/O. Each rank directly
  writes its disjoint extent; payloads never travel through rank 0.
- Empty aligned input chunks publish an unchanged prefix from the reader.
- Close collectively checks the number of published chunks. Failed operations
  return errors to the application's existing fatal-error path.

The prefix slots consume eight bytes per nominal chunk on rank 0. The writer
keeps only its current packed chunk, using the existing bounded pipeline/SAM
storage. There is no unbounded completed-SAM reorder queue. Offset arithmetic
is checked against the signed 64-bit file-offset limit.

Monotonic tickets are essential: each rank's processing/output FIFO is ordered,
so the earliest unpublished chunk cannot be waiting behind a later chunk on
that same rank. Arbitrary distributed queue order would invalidate this
argument. An explicit `SWBWA_MPI_TICKET_MODE=distributed` is rejected.

This first implementation uses prefix polling with a scheduler yield, not
an artificial delay or a per-chunk collective barrier. It may create substantial
RMA traffic and suffer head-of-line blocking on a heavy chunk. This is a real
limitation to profile, not a recommended way to slow another configuration.
Future batched prefix publication can reduce polling without changing output
semantics. `ordered prefix polls` and `ordered prefix wait and publication`
are reported with output diagnostics.

## Validation

Run `python3 tests/test_output_ordered.py` using a native MPI installation.
The test uses three ranks, deliberately delays an early chunk, checks actual
file bytes in input order, includes empty files and empty aligned chunks,
tests exact/fast indexing and tail refinement, and repeats open/close to check
truncation. It covers one/six packing workers and `OUTPUT_RMA_ONLY=1`.

`OUTPUT_RMA_ONLY=1` preserves prefix ordering, all SAM generation/packing, and
RMA, but does not create a file. Extent offsets must be contiguous **in chunk
ID order**, not merely nonoverlapping. Its sampled SAM fingerprint is not a
full-file checksum or a substitute for testing the real writer.

The generic streaming output API has no chunk ID and is rejected in ordered
mode. Use `swbwa_output_write_chunk()`. No SAM header is emitted in this mode.
Native tests alone do not establish Sunway performance or target ISA correctness.

For actual alignment validation, compare a small real FASTQ subset against a
non-MPI sequential reference using `tests/check_ordered_sam.py`. Enable
`MPI_EXACT_READ_INDEX=1`, retain real writes (`OUTPUT_RMA_ONLY=0`) and fix the
insert distribution for PE. The checker compares every record byte without
sorting; run both PE/SE and synchronous/pipelined variants. Fast read indexing
may change equal-score tie breaking, so an exact-reference MD5 check must not
silently compare exact and fast modes. Keep these correctness runs separate
from performance runs that omit pwrite.
