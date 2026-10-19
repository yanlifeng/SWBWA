# Parallel Gzip FASTQ Input

The public C interface is `include/swbwa_gzip.h`, implemented by `input.cpp`.
It preserves input order, emits matching record counts for paired inputs,
and never reads or changes alignment's `n_processed`. Returned allocations
are writable, NUL-terminated, and independently owned by the caller (`free`).
Byte lengths include the final LF and exclude the sentinel. The caller must
drain to EOF to validate the complete stream; errors are sticky.

Input is gzip containing unwrapped, four-line, LF-terminated FASTQ. Missing
final LF is deliberately rejected, not repaired: newline counts alone cannot
distinguish a complete unterminated quality string from a truncated one.
Incomplete groups of four lines, truncated/corrupt gzip, unequal mate record
counts, and nonzero trailing garbage are errors. Concatenated members, empty
members, and final zero padding are supported. FASTQ headers, sequence/quality
lengths, and matching mate names are left to the existing formatter; this is
not a general FASTQ validator or a multiline FASTQ parser.

## Integration

```c
swbwa_gzip_reader *reader = swbwa_gzip_open(path1, path2, chunk, workers, backend);
/* On NULL: swbwa_gzip_error(NULL). Backend is "rapid" or "zlib". */
/* Handle nonzero status as an input setup error. */
swbwa_gzip_set_batch_bytes(reader, fastq_batch_bytes);
char *r1, *r2;
int64_t n1, n2;
int status;
while ((status = swbwa_gzip_read_next(reader, &r1, &n1, &r2, &n2)) == 1) {
    /* Required before using the legacy shared-offset PE formatter. */
    if (swbwa_gzip_check_pe_layout(reader, r1, n1, r2, n2, cpe_count) != 0) {
        free(r1); free(r2);
        status = -1;
        break;
    }
    /* Hand owned buffers to the existing alignment pipeline. */
}
/* status == -1: swbwa_gzip_error(reader); status == 0: verified EOF. */
swbwa_gzip_report(reader, stderr);
swbwa_gzip_close(reader);
```

SE uses `path2 == NULL`; all four `read_next` output arguments are still
required. `swbwa_gzip_next` is an alias. Error/EOF output buffers and byte
counts are NULL/zero. One reader is sequential and is not thread-safe.
The API does not read `SWBWA_GZIP_READER`; the application passes its value.

`chunk_bytes` is the compressed discovery interval (typically 256 KiB).
`swbwa_gzip_set_batch_bytes` independently sets the decoded FASTQ batch target
before the first read (default 16 MiB of R1). The application should pass
its existing per-CG batch size times CG count. A complete record longer than
the target is returned alone; otherwise R1 batches end at/before the target.
PE refills R2 until the R1-selected record count is available; R2 byte lengths
do not shrink R1 batches and may exceed the R1 byte target.
A record can cross multiple refills. Memory is bounded by codec batches, queued decoded views,
and the longest unfinished record, not by the complete input file. Codec
workers are 1..384; chunks are 64 KiB..64 MiB. The existing engine tuning
variables are `SWBWA_GZIP_PREFETCH_MIB` (0..256, default 64) and
`SWBWA_GZIP_SYMBOL_BUDGET_MIB` (512..2048, default 2048).

## Pair Layout

`read_next` permits different byte lengths/layouts and emits equal record
counts; that guarantee alone does not make a batch safe for the existing
shared-offset CPE formatter. Call `swbwa_gzip_check_pe_layout` immediately
after reading and before modifying/submitting buffers, with the formatter's
actual CPE count, not codec workers. It rejects unequal batch byte lengths
and checks equal record ranks at every R1-selected CPE boundary. It mirrors
the existing `ceil(bytes / CPEs)` division, 1024-byte lookbehind, and
`get_next_fastq` boundary selection. Below `CPEs * 1024` bytes the formatter
uses a single partition. Interior record lengths need not match where the
two independent parsers consume the same complete record range.

Page indexes locate the boundary lines without walking every record. The
rapid resolve kernel counts newlines while packing bytes; the host uses
counts for full pages, scans partial pages/boundaries only, and makes one
copy of each delivered byte into caller-owned memory. Pending views are
not rescanned or copied on each refill. Serial zlib/fallback has no decoded
page counts and constructs its index with one newline pass per decoded
buffer; there is no additional full FASTQ parse in the wrapper. The
metrics distinguish decoder work and returned-buffer copy costs.

## Build Boundary

The application Makefile selects this engine with `GZIP_INPUT=1`; `build.sh`
retains two-pass cross-segment layout discovery and relocation. Runtime
ownership is implemented in `src/host/bwamem.c`; the CPE entry lives in
`src/slave/slave.c`. Set `GZIP_INPUT=0` for an application without codec objects.
Host objects: `input.cpp`, `gzip_stream.cpp`, `rapid_reader.cpp` (C++11).
Codec objects: `codec_job.c`, `rapid_firstpass.cpp` (C++17), plus exactly one
backend. Portable tests link `backend_host.c`; Sunway builds additionally
compile `backend_swbwa.c` as host code, providing `swbwa_codec_launch/init/finish`.
It calls `swbwa_codec_runtime_launch` from `include/swbwa_codec_runtime.h`.
The application's runtime dispatches `backend_cpe.c`'s `swbwa_codec_worker`
through its own runtime/crossentry.

Use `SWBWA_CODEC_RAPID=1` for rapid support. The CPE first pass supports
`SWBWA_CODEC_CPE=1`, `SWBWA_CODEC_LDM=2`, retaining the existing allocation
failure/stack-headroom fallback in `rapid_state.h`. The application owns
init/finish and serialization of codec launches against complete alignment
operations. Reader open/close never initialize or halt the shared runtime.
The alignment formatter's existing maximum reads per batch still applies;
very short records may require reducing `-K` even when the byte budget fits.

CPE objects may use `Z_PREFIX` and `vendor/zlib` with `adler32.c`, `crc32.c`,
`deflate.c`, `trees.c`, and `zutil.c`. Host decoding/fallback links system
zlib and pthreads. The inherited `SwbwaGzipWriter` dependency remains inside
the engine; this API adds no writer or SAM gzip output workflow.

Run standalone host tests without the application Makefile:

```sh
python3 tests/test_gzip_input.py
python3 tests/test_gzip_input.py --zlib-only
python3 tests/test_gzip_input.py --sanitize
```

## Provenance

The codec engine and minimal vendor subset were adapted from RabbitBio/SWQC
(`https://github.com/RabbitBio/SWQC`), `codec/` at revision
`0d1bcf130010ba8d243a787dc342e27458a9361c`. The SWQC MIT license is preserved
in `LICENSE-SWQC`; upstream rapidgzip MIT/Apache and zlib license texts and
source notices remain under `vendor/`. `vendor/NOTICE` records the upstream
revisions and SWQC's Sunway-specific adaptations.

SWBWA changes: symbol/include-guard/environment namespace conversion,
runtime-independent CPE entry, strict serial gzip member validation, indexed
line-rank queries, and the owned-buffer C FASTQ/pair-layout API and tests.
No claim is made that these adaptations are upstream rapidgzip behavior.

## MPI input

`include/swbwa_mpi_gzip.h` / `mpi_input.cpp` connect the distributed decoder
to the existing three-stage alignment pipeline. Build with `USE_MPI=1`,
`GZIP_INPUT=1`. `EXEC_MODE=single` gives each rank 64 CPEs and allows six
ranks per node. `cgs` and `cgs_cross` give each rank a full-chip runtime;
do not launch multiple full-chip ranks on the same node. Standard modes use
the existing spawn/join and MPI-progress completion path; `cgs_cross` uses
the relocated dispatcher. Both serialize codec launches against alignment.
In a threaded pipeline, standard spawn/join requests from the reader are
executed by the processor/main thread while it waits for input. The input
queue wakes that thread without polling; its mutex is released during the
launch. This preserves athread ownership while file reads, framing and SAM
writes can still overlap alignment. Codec wait time includes this handoff,
not just the runtime mutex. The cross dispatcher retains its direct path.
MPI thread support
must satisfy the existing SWBWA `MPI_THREAD_MULTIPLE` requirement.

`mpi_decoder.cpp` uses bounded compressed windows and multiple local CPE jobs
per collective epoch. Job metadata establishes a contiguous decoded chain;
32 KiB dictionaries propagate in stream order and CPEs resolve history
references. Ranks retain their local decoded bytes. CRC/ISIZE are validated
for every gzip member, including concatenated members. No preliminary full
decompression or external index is required. This does not implement arbitrary
independent seeks into a gzip stream.

MPI gzip input assumes files remain unchanged throughout the run. The decoder
uses one initial `fstat` to obtain the file size and require a regular file;
it does not compare timestamps, monitor file identity, or restat between
epochs. Stale filesystem timestamps therefore cannot trigger input-change
errors. Collective open still checks the path, size, gzip header/trailer and
decoder configuration. Read errors, truncation, malformed gzip data and member
CRC/ISIZE failures remain fatal. These checks do not guarantee detection of
concurrent file modifications; callers must keep the input immutable.

`mpi_input.cpp` uses exact global newline prefixes and the decoder's sparse
4 KiB newline counts, not an `@` search or a full MPE FASTQ scan. Only partial
pages and requested record boundaries are scanned. At most one partial record
follows each decoded boundary, capped at 1 MiB. Complete records travel as bulk
bytes plus sparse newline metadata, not as individually serialized records.
Local destinations bypass MPI packing and self-exchange; large copies use the
existing MPE helper pool without waiting behind another pool user.

Final output chunks are independent, caller-owned mutable allocations. A
`swbwa_fastq_partition_t` array accompanies them through Stage 2. R1 boundaries
are selected by byte position and rounded to complete record ordinals; R2
tries the same byte position, then uses the newline index to correct it when
necessary. Equal byte sizes or matching boundary names alone are not proof of
pair alignment. The CPE formatter parses each partition once, validates every
four-line record and compares every pair's exact name tokens. It reuses known
partition counts instead of counting all newlines again. The old shared-byte-
offset guard is only needed by the non-MPI reader, not this path. Free the
partition array after formatting, alongside the ordinary batch lifecycle.

All ranks call `next` in the same wave order. A READY result may contain zero
records on a rank; the pipeline must forward that empty token without launching
alignment or writing an output chunk. EOF is returned collectively only after
both mate trailers and final record counts have been checked. On a local parse,
allocation or MPI failure the adapter aborts the communicator, rather than
leaving peers blocked in a collective. Normal close is collective after the
pipeline has drained.

The first decoded R1 records estimate a fixed records-per-chunk target from
`-K`; chunk `i` belongs to rank `i % ranks`. Read indices are always exact and
computed online, independent of `MPI_EXACT_READ_INDEX`. The plain FASTQ static/
dynamic scheduler and tail refinements are bypassed. Supported output modes
are `split`, `single_unordered`, `single_ordered` and `discard`.
For `single_ordered`, the output layer holds an online `(next chunk ID, next
byte offset)` prefix in an RMA window instead of requiring a precomputed chunk
count. The unique owner of the next ID publishes its byte prefix before its
successor ID. Empty READY tokens in the final wave do not publish chunks.
After validated EOF and writer completion, every rank supplies the exact global
chunk count; close validates both the final count and reserved byte total.
This does not introduce dynamic gzip input scheduling. Static ordered builds
are supported for gzip only; plain FASTQ still requires dynamic global tickets.
The report separates decode, framing/boundary exchange, record redistribution
and final assembly/partition planning. `reader_local` subtracts adapter MPI
calls and READY/EOF barriers from local work; decoder reports separately expose
pread, local work, CPE execution, runtime-queue wait and MPI phases. Timers are
accumulated wall time, not CPU time summed over helper threads. MPI time includes
waiting for peers (which may still be aligning). In has1 mode there is no local
processor overlap; still exclude decoder runtime wait when comparing Reader
work. Full FASTQ validation remains included in Stage 2 formatting, so report
that subphase too rather than claiming moved work as a total speedup.

`OUTPUT_RMA_ONLY=1` retains SAM generation, packing and ordered RMA publication
but does not create a SAM file. Its default validation samples 100 reads per
chunk. `SWBWA_OUTPUT_SAMPLE_READS=0` instead hashes every packed SAM blob,
counts all SAM lines and reports a full packed-byte CRC32 for each chunk.
Combining these CRCs in chunk order (with their byte lengths) checks the entire
ordered byte stream against the CRC reported by non-MPI streaming MD5. This
is not a distributed MD5 and is not a collision-free proof of equality.
Full validation cost is included in Stage 3; report the sampling policy with
performance measurements. Normal builds keep physical output enabled.

The record target is an estimate, not a hard byte limit for every later chunk.
Substantially longer records later in a file can exceed the formatter capacity
and cause an explicit error; chunks are not adaptively split in this first MPI
implementation. The reported queue limit is per mate, not total process RSS:
decoder windows, record metadata, exchange buffers and in-flight pipeline
batches also occupy memory. MPI redistribution can therefore cost more memory
and host work than the single-process gzip path.

Host tests exercise the real decoder with a host codec launcher; they do not
replace native CPE/runtime validation:

```bash
python3 tests/test_mpi_gzip_decoder.py --work /tmp/swbwa-mpi-decoder-test
python3 tests/test_mpi_gzip_input.py --work /tmp/swbwa-mpi-input-test
python3 tests/test_mpi_gzip_pipeline_contract.py
```
