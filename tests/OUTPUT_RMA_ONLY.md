# Output RMA-only benchmark

Build with `USE_MPI=1 MPI_INPUT_MODE=dynamic OUTPUT_MODE=single_unordered
OUTPUT_RMA_ONLY=1`. The default `OUTPUT_RMA_ONLY=0` leaves production behavior
unchanged. Both compilation and Make reject unsupported output/input modes.

This benchmark generates all alignments and SAM strings. It keeps the normal
output packing and reservation work, but no SAM file is created, truncated or
written. The alignment pipeline now submits complete input chunks: each
nonempty chunk reserves its full SAM byte length with one existing MPI
fetch-and-add/flush pair. The packing buffer grows beyond its initial size
when necessary. Legacy streaming API callers still use the buffered flush
schedule. Input RMA and output RMA therefore coexist, without disk
write traffic. This is not a filesystem throughput benchmark or a claim about
hardware RDMA offload.

The first 100 input reads of every chunk are sampled (50 pairs for PE). Each
sample hashes the complete SAM string for that read, including supplementary
records. End-of-run logs, printed in rank order, include chunk ID, byte range,
read count, full SAM byte count, sample count/bytes and sum/XOR hashes. A final
per-rank summary includes submitted/reserved bytes and reservation counts.
Require unique complete chunk coverage and matching ranges/counts/samples in
cross-run checks. This is sampled consistency checking, not full SAM accuracy
validation. A matching total byte count cannot detect all unsampled changes.

With exact read indices disabled, dynamic chunks use their aligned input byte
offset as the deterministic read-index seed. Comparison requires the same
actual chunk boundaries, read order and alignment parameters. It does not imply
equality with standard BWA's record-index seed or with different chunk sizes.

Use `-v 4` for stage, scheduler and output timing details. All samples are
collected by the existing stage-3 writer; no additional communication occurs
per chunk. Ordered reporting at shutdown is outside the pipeline measurement.
The pipeline timer now stops when all pipeline threads and their final output
flush finish; it excludes scheduler/output diagnostic printing and MPI window
teardown. Older logs timed those operations as part of the pipeline, so do not
mix their total-time column into this experiment's scaling curves.
For launchers that interleave process output at byte boundaries, set
`SWBWA_RANK_REPORT_PREFIX=/path/to/run`. After the pipeline timer stops, each
rank temporarily redirects its final reports to `run.rank000000.log`, etc.,
then restores stderr before normal cleanup. The original launcher pipe stays
open throughout. The directory must already exist, and existing rank files
are not overwritten. Keep the launcher log as well: initialization,
configuration and the final command line remain there. This option is only
active in the RMA-only benchmark and does not write SAM records.
Run `python3 tests/test_output_rma.py` for the standalone real-MPI regression.
See [CHUNK_OUTPUT.md](CHUNK_OUTPUT.md) for the chunk contract and optional
`HOST_MPE_THREADS=6` helper team. Older 64 MiB streaming benchmark logs are a
different output-granularity baseline and must be labelled accordingly.
