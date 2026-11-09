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
zlib and pthreads. SAM output uses the separate `output.c` scatter/gather engine,
not the inherited staging writer in `gzip_stream.cpp`.

## Parallel SAM Compression

`GZIP_OUTPUT=1` (the Makefile default) enables `.gz` output for split,
single-unordered and single-ordered modes. Input compression is independent:
`GZIP_INPUT=0 GZIP_OUTPUT=1` is supported. Other output names remain uncompressed.

Each compression task produces a complete RFC 1952 gzip member with its own CRC32
and ISIZE. Concatenating members produces an ordinary gzip file, not BGZF or a
private container. Member boundaries may split a SAM line; decompression restores
the identical byte stream. No history window is exchanged between ranks.

The engine borrows SAM addresses and lengths until synchronous compression returns.
MPE helpers measure lengths and build prefix metadata; CPEs consume those spans
directly, without an MPE-sized SAM concatenation. Compressed members occupy reused
slab slices and are sent through writev/pwritev without a second packing copy.
Each CPE reuses zlib state across its jobs within one launch. A bounded 32 KiB LDM
allocation holds small zlib state/hash/pending arrays when available; larger arrays
use a reusable, 64-byte-aligned 160 KiB main-memory arena per CPE. Only exhaustion
falls back to libc allocation. Alignment and codecs share runtime
ownership, and native spawn/join remains on the processor/main thread.

By default (`SWBWA_GZIP_PREPARE_OUTPUT=1`), the processor compresses each batch
immediately after alignment, before handing it to the writer. The writer then owns a
prepared compressed batch until its positioned write completes. Encoding and
physical writes can overlap without reusing a slab that is still being written.
There are `SWBWA_PIPELINE_BUFFER_COUNT` independent result slots, currently five;
each retains its high-water capacity. This bounds outstanding batches, not their
total bytes, so larger `-K` values also increase potential compression-buffer
memory. Encoder scratch storage is shared, not replicated per slot.

Preparation is reported separately as `gzip preparation between stages 2 and 3`.
It is excluded from the alignment Stage 2 timer and is not eliminated work.
Compare complete Pipeline time as well as this timer when evaluating the option.
The writer-requested compression path remains available with the option set to 0.

The default level is 1. `SWBWA_GZIP_LEVEL=1..9` changes it. Automatic member sizes
target two tasks per configured CPE, bounded to 16..256 KiB; there is no 64-task
batch cap. `SWBWA_GZIP_BLOCK_KIB=1..4096` selects a fixed size, and `0` restores
automatic sizing. Small inputs need not occupy every CPE. Single-CG builds use up
to 64 workers; cgs/cgs_cross use up to 384. Six-MPE builds use their existing bound
helper pool for descriptor preparation, falling back to the caller when busy.

The MPI output layer reserves one compressed extent per input chunk. Unordered
output reserves extents in completion order; ordered output publishes prefixes
by input chunk ID. Actual pwritev calls may complete in any order.
Six-MPE builds can partition compressed vectors over a separate bound I/O pool for
disjoint positioned writes with `SWBWA_GZIP_PARALLEL_WRITE=1` (default off).
The I/O pool uses the same allocated MPEs but has its own dispatch lock, so a
blocked filesystem write does not reserve the compute helpers needed by Reader.
Per-thread syscall time is a
worker sum; `[SWBWA gzip sink ...] write_wall_s` is elapsed write time.
Ordered-prefix polling yields after eight unsuccessful probes, using an absolute
monotonic deadline resilient to interrupted sleeps. `SWBWA_ORDERED_POLL_US=0..1000`
sets the optional pause (default 0, yield without sleeping). Actual wakeup latency depends
on the system timer resolution. No CPE/helper mutex is held during the pause.
`SWBWA_GZIP_PREFIX_AHEAD=1` is an opt-in ordered-gzip-output experiment. After Processor
compresses a batch, it submits only its ID and compressed length to a bounded
metadata worker. The worker publishes the ordered prefix independently of
Writer's previous disk write; Writer still waits for its own offset before
writing its existing compressed vectors. There is no additional SAM copy or
payload queue. The worker shares the last allocated MPE in a multi-MPE build
and requires MPI_THREAD_MULTIPLE. An announced batch must be written, not
released without writing. Both plain and gzip inputs can use it; unordered,
split and unprepared output keep their original path. Each rank must submit
strictly increasing chunk IDs, as required by the ordered Writer as well.
While idle, the metadata worker relies on the existing MPI progress thread;
only when that thread is inactive does it periodically enter MPI_Iprobe itself.
These fallback probes and wall time are reported separately. No message is
consumed and no CPE lock is held.
With prefix-ahead enabled, the ordinary reservation timer belongs to this
background worker and overlaps Processor/Writer. Use its `writer_wait_s` to
measure the Writer's exposed prefix stall; do not add the background timer to
Stage 3 or Pipeline elapsed time.
An empty output is still a valid empty gzip file. Split names are
`name.rank000000.sam.gz`.
Discard and streaming-MD5 modes reject `.gz` output instead of silently writing
another format. RMA-only profiling still runs compression, reserves compressed
lengths, and fingerprints the original SAM bytes, but does not create a file.

`-v 4` file-output counters and the per-rank `[SWBWA gzip output ...]` line separate
MPE preparation, CPE dispatch/execution, runtime wait, compressed byte volume and
physical writes. Codec runtime wait is not compression work; a pipelined run may
include alignment/runtime contention in stage times.

Local regression: `python3 tests/test_gzip_output.py -v`. It covers all output
modes, concatenated member validation, empty/idle ranks, oversized SAM records,
reopen, compression settings, partial writes and EINTR. Set
`SWBWA_TEST_SANITIZER=address,undefined` to instrument the production C paths.

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
By default, in a threaded pipeline, codec requests from the reader and writer are
executed by the processor/main thread at batch boundaries, including queue
waits and final draining. Requests wake that thread without polling; the queue
mutex is released during the launch. This preserves athread ownership while file reads, framing and SAM
writes can still overlap alignment. Codec wait time includes this handoff,
not just the runtime mutex. Cross-mode codec requests use the same ownership
handoff before entering the relocated dispatcher.
`SWBWA_GZIP_CPE_DISPATCH=reader` is an experimental cross-mode alternative:
initial native runtime setup still runs on the main thread, then the decoding
caller publishes directly to the persistent dispatcher under the runtime mutex.
Alignment retains its whole-batch lock and compression stays on the processor.
Standard modes keep main-thread athread launches even when `reader` is requested;
the startup log reports both requested and effective policies. The default is
`processor`. At verbosity 4, per-operation reports separate handoff, mutex wait,
and runtime execution; direct dispatch is not evidence of concurrent CPE kernels.
`SWBWA_GZIP_CPE_DISPATCH=hybrid` additionally tries the existing MPE helper pool
when decoding finds the CPE runtime occupied. It executes the same discover and
resolve implementations on independent job buffers, without touching CPE state.
A busy helper pool falls back to the regular CPE path. This is experimental and
applies to standard mode too (one MPE per rank); it does not reserve extra cores.
Compression still uses CPEs on the processor. Reports separate `run_s` (CPE)
from `mpe_s` (helper execution); decoder `discover_run_s` / `resolve_run_s` include
whichever backend ran. MPE decode may overlap alignment, so do not add these
durations to Pipeline time. The default remains `processor`; the experimental
hybrid policy is opt-in, and reducing codec waits alone does not guarantee a
faster pipeline when input or output is the bottleneck.
MPI thread support
must satisfy the existing SWBWA `MPI_THREAD_MULTIPLE` requirement.

`mpi_decoder.cpp` uses bounded compressed windows and multiple local CPE jobs
per collective epoch. Job metadata establishes a contiguous decoded chain;
32 KiB dictionaries propagate in stream order and CPEs resolve history
references. Ranks retain their local decoded bytes. CRC/ISIZE are validated
for every gzip member, including concatenated members. No preliminary full
decompression or external index is required. This does not implement arbitrary
independent seeks into a gzip stream.

A syntactically valid speculative candidate need not join the proven prefix.
Before abandoning the rest of an epoch, the decoder can repair the job at that
gap, using the exact start bit established by the preceding job. The repair
starts with history markers, not an invented dictionary; the existing ordered
dictionary propagation and RESOLVE provide the actual bytes. No extra input
read or symbol-buffer allocation is required. `SWBWA_GZIP_GAP_REPAIR_LIMIT`
selects 0..8 repairs per epoch (default 8, 0 disables). All ranks must agree.
If a repair needs more output capacity, only the old proven prefix is committed;
normal next-epoch anchor growth still handles long blocks. CRC/ISIZE checks are
unchanged. `repair.total_s` is inclusive of backend/MPI time, not an extra phase
to add to `next_s`; overwritten candidates count as discarded jobs.

Compressed bulk reads use a persistent, separately synchronized I/O pool.
`SWBWA_GZIP_INPUT_READERS` selects 1..`SWBWA_HOST_MPE_THREADS` readers (default:
all allocated MPEs). Regions of at least 1 MiB use disjoint `pread` ranges;
small header/trailer reads remain serial. These threads never call MPI or CPEs,
and cannot hold the compute-helper lock during filesystem waits. Exact overlap
between successive rank-local compressed windows is copied into the bounded
prefetch allocation; `SWBWA_GZIP_REUSE_INPUT=0` disables that reuse. It never
reuses speculative decoded bytes. Reports separate source-request bytes from
cache-reused bytes; neither is a measurement of physical disk traffic.

For threaded pipelines, `SWBWA_INPUT_QUEUE_CHUNKS=1..8` overrides the bounded
input queue (default 1 for MPI, 4 without MPI). It is independent of the gzip
ready-window size below. Extra queued chunks can overlap refill with alignment,
but also increase resident memory and advance task claims; deeper is not always
faster. Writer/SAM ring capacities and buffer lifetime rules do not change.

Decoder concurrency and budgets follow the CPE count in both static and
dynamic modes. A 64-CPE rank allows 64 jobs with a 512 MiB decoder-buffer budget
and 128 MiB decoded-output limit per epoch. A 384-CPE rank allows 384 jobs with
3 GiB and 768 MiB limits, respectively. The compressed stride stays 256 KiB and
the initial symbol capacity stays 2 MiB per job. Limits are demand-driven, not
preallocations; short inputs, rejected candidates and high-compression capacity
retries can still reduce useful parallelism. Reports expose `initial_batch_jobs`
and the actual maximum DISCOVER/RESOLVE batch widths so the configured CPE count
is not mistaken for the number of jobs doing useful work.
When fewer compressed strides remain than the full job grid, each rank's
stripe shrinks to share the remaining strides across ranks. This avoids giving
all of a small file or tail to the first few ranks. It cannot create additional
independent tasks when too few DEFLATE blocks remain.

The per-mate FASTQ queue budget is the larger of four decoded-output limits and
three alignment byte targets (512 MiB or 3 GiB at the default targets). It is
separate from decoder buffers, transient exchanges, ready-window payloads and
in-flight alignment allocations; these figures are not process RSS limits.
Increasing decode concurrency does not change the `-K` alignment target or the
number of dynamic ready chunks.
Compressed windows are filled directly by `pread`; only the bounded decoder
lookahead guard is zeroed, avoiding an extra host write of the entire window.

After a validated epoch, each rank asynchronously reads the exact compressed
stripe required by its next epoch. The I/O worker only calls `pread`; it never
enters MPI or launches CPE work. Its allocation is transferred directly to the
decoder. Unlike sequential file readahead, this does not intentionally fetch
stripes belonging to other ranks. The existing `SWBWA_GZIP_PREFETCH_MIB` setting
limits background bytes per request (MPI default 128 MiB, range 0..256); any
remaining bytes are read on demand, and 0 disables prefetch. No full-file scan
or startup body read is added. Retained output and the prefetched compressed
allocation must fit the decoder budget. Full-chip builds bind the mate I/O
workers to auxiliary MPEs and verify their placement.

`foreground_input_s` (also the legacy `io_s`) measures reader-visible input
time, including waiting for a prefetch. `source_io_s` measures actual `pread`
time, with `prefetch_io_s` its background subset. Background time overlaps other
work and must not be added to Reader time. Source byte/call counts include the
consumed prefetch requests; validated EOF leaves no pending request.

Each MPI decoder Source reuses one lazily started I/O thread, with at most one
pending request. It does not create a new thread for each compressed window.
Worker exceptions are delivered through the pending future, and close joins the
worker before releasing the file descriptor. The worker never calls MPI or uses
the CPE runtime; native MPE placement is retained for its lifetime.

The optional MPE decoder compiles with its distance/history buffers embedded
(`SWBWA_CODEC_MPE=1`, `SWBWA_CODEC_LDM=0`). Its state uses malloc-backed manual
alignment and frees the original pointer. This avoids the native aligned
new/delete path retaining ordinary MPE memory under xmalloc. CPE scratch and
its LDM placement remain unchanged.

DISCOVER also supplies a byte dictionary when its last 32 KiB contain no
unresolved references, using the same codec fast path as the non-MPI reader.
The MPE reuses that dictionary instead of rescanning the 16-bit symbols.
Unknown or short tails still propagate through the checked history-repair
path. `known_tail_jobs` and `repaired_tail_jobs` count these accepted-job paths.

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
Remote record descriptors retain their source buffers until exchange completes.
Their headers, sparse newline metadata and payloads are written directly into
the final MPI send allocation, avoiding intermediate per-destination payload
buffers and an additional flattening copy. Send and receive buffers are not
zero-filled when every transmitted or consumed byte is explicitly initialized.

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

In static mode, all ranks call `next` in the same wave order. A READY result may contain zero
records on a rank; the pipeline must forward that empty token without launching
alignment or writing an output chunk. EOF is returned collectively only after
both mate trailers and final record counts have been checked. On a local parse,
allocation or MPI failure the adapter aborts the communicator, rather than
leaving peers blocked in a collective. Normal close is collective after the
pipeline has drained.

The first decoded R1 records estimate a fixed records-per-chunk target from
`-K`; decoding initially stages chunk `i` at rank `i % ranks`. Read indices are always exact and
computed online, independent of `MPI_EXACT_READ_INDEX`. The plain FASTQ static/
dynamic scheduler and tail refinements are bypassed. Supported output modes
are `split`, `single_unordered`, `single_ordered` and `discard`.
For `single_ordered`, the output layer holds an online `(next chunk ID, next
byte offset)` prefix in an RMA window instead of requiring a precomputed chunk
count. The unique owner of the next ID publishes its byte prefix before its
successor ID. Empty READY tokens in the final wave do not publish chunks.
After validated EOF and writer completion, every rank supplies the exact global
chunk count; close validates both the final count and reserved byte total.
Static ordered builds are supported for gzip only; plain FASTQ still requires
dynamic global tickets. Dynamic gzip scheduling is described below.
`SWBWA_GZIP_DECODE_DRAIN=1` restores the optional full alignment-queue drain
before the next dynamic decode epoch; the default is 0. Disabling that early
drain does not remove CPE runtime ownership, READY-window completion, collective
epoch synchronization or EOF draining. It permits overlap but cannot make
alignment and decoding use the same CPE runtime simultaneously.
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

### Dynamic Decoded-Chunk Scheduling

`MPI_INPUT_MODE=dynamic` retains collective decompression but decouples alignment
ownership from the decoding rank. `SWBWA_GZIP_READY_CHUNKS` selects 1..16 staged
chunks per rank (default 4); all ranks must use the same value. A refill stages
up to `ranks * ready_chunks` complete chunks, including independent mate buffers,
exact global read ordinals and CPE partition metadata. No whole-file scan or
compressed-stream random-access index is introduced.

A refill reserves its first wave of chunks, one per rank, for their staging
owners. These allocations transfer directly to the caller, without another
pack or copy. A rank-0 64-bit `MPI_Fetch_and_op`/`MPI_Win_flush` ticket assigns
the remaining increasing global IDs. A rank must consume its reserved chunk
before claiming a ticket: all reserved IDs precede the ticketed suffix, so
local output order remains increasing. With `SWBWA_GZIP_READY_CHUNKS=1`, every
chunk is reserved locally; use a larger window for dynamic load balancing.

`SWBWA_GZIP_EARLY_FIRST=1` is an opt-in delivery policy: return the owner's
reserved first chunk as soon as it is assembled, before staging the rest of
the window. The next call prepares and publishes the ticketed suffix using
the same chunk IDs and RMA protocol. This lets alignment of the first chunk
overlap Reader preparation without increasing the alignment input queue.
All ranks must select the same policy. The default is 0 for comparison.

Only ticketed data is packed into an ordinary MPI window with relative byte
displacements. Remote claims use `MPI_Get`/flush, while local ticketed claims
copy from that window into caller-owned buffers. The window and allocation
are reused while capacity suffices, and rebuilt collectively only on growth.
Capacity is rounded up to MiB boundaries to absorb small record-length changes.
Sources are immutable until all Gets are flushed and every rank reaches the
refill-entry barrier. `MPI_Win_sync` publishes the next generation before
the reset ticket is made available. Returned buffers remain owned by the caller
even after subsequent refills. PE buffers are always claimed together. This deliberately
avoids dynamic-window absolute addressing: a standalone native SWCH test passed
ordinary-window Get but failed dynamic-window Get with `SRC_ADDR_ERR`.

With one MPI rank there is no ownership imbalance to schedule. The reader
directly returns decoded chunks, bypassing scheduler RMA, packing and staging
windows. Codec/runtime locking and output-layer MPI progress remain unchanged.

`next()` is not collective per chunk in dynamic mode. It returns only nonempty
READY batches; ranks may receive different numbers of them. It is unsafe for a
caller to add an Allreduce/Barrier after every READY. Global EOF is known only
after trailer validation and depletion of the final ready window. Chunk IDs
increase on each rank, which preserves the ordered writer's FIFO contract.

After Stage 3 consumes each chunk, `swbwa_mpi_gzip_complete(reader, id)` releases
its local in-flight credit. This does not free FASTQ/SAM allocations or call MPI.
Before an actual decoder advance, the reader waits for all outstanding local
credits only when `SWBWA_GZIP_DECODE_DRAIN=1`. This optional drain includes
chunks delivered by earlier ready windows. The default (`0`) allows decoding
to overlap prior alignment/output. Processor-side codec requests are serviced
while waiting on either input or output queue space, and during final writer
draining; a full output queue must not prevent a pending codec request from
running. The immutable source-window lifetime is protected independently by
flushed Gets and the refill barrier. `has1` also calls complete, normally before
its next reader call. Static and single-rank
completion are no-ops. EOF can now arrive with output still in flight; join the
pipeline and finish all completion callbacks before closing the reader. Close
rejects outstanding distributed credits. The existing MPI progress thread stays
active throughout the pipeline.

The standard-mode reader requests codec launches on the processor thread using
pthread-managed thread-specific storage. This preserves its queue association
across native alignment launches; the earlier compiler-TLS pointer was observed
to change to an invalid queue address after the first alignment window.

`SWBWA_GZIP_PREPARE_AHEAD=1` optionally gives distributed dynamic input a
bounded producer thread in the no1 pipeline. It prepares at most one READY
window ahead (plus one in-flight chunk) while the consumer claims the previous
window. The producer owns decoder/framing state, the consumer owns RMA windows
and tickets, and their collectives use separate communicators. The producer
inherits the pipeline's processor-handoff context, so native CPE launches keep
their owner-thread requirement. EOF joins the producer before destroying the
pipeline context. The default is 0; single-rank, static and has1 execution keep
the synchronous path. All ranks must agree on the option.

The `prepare_ahead` report separates producer wall time, producer backpressure,
consumer preparation waits and peak queued FASTQ bytes (not total RSS). With
this option, background decode/framing/layout timers overlap `next_s` and are
not children of it or `refill_s`; do not sum them into foreground Stage 1.
Midstream reports include only consumer-owned scheduler fields. Full decoder
statistics are reported after the producer has joined at EOF.

`SWBWA_GZIP_PARALLEL_HISTORY=0|1|2` selects the dictionary propagation strategy
(default **0**, identical on every rank):

| Value | Strategy | History propagation critical path | Summary storage/rank |
| --- | --- | --- | --- |
| 0 | Serial owner preparation and broadcasts | `O(P J B)` | `O(B)` |
| 1 | Local summary, Allgather, local linear prefix | `O(J B + P B)` | `O(P B)` |
| 2 | Local summary, ordered binary tree up/down | `O(J B + B log P)` | `O(B log P)` |

Here `P` is the rank count, `J` the maximum accepted jobs/rank, and `B=32768`.
These are history-work bounds, not whole-decoder speedup promises. Actual local
summarization traces only necessary symbols rather than always touching `J B`.

Modes 1 and 2 use the same local summary construction.
Each rank composes its accepted job tails into one 32-KiB symbolic dictionary:
each 16-bit element is a literal or a reference to that rank's incoming window.
Only the final window and the earlier symbols it references are visited; a
marker always jumps before its job, so tracing terminates. This avoids composing
32 KiB for every job just to produce the rank's final summary. A bounded,
direct-mapped cache reuses resolved positions in long shared reference chains;
replacement can cause a cache miss but cannot change the result. Its storage
is charged to the decoder memory budget. Regression tests
compare the result against full sequential symbolic composition, including
chains of tiny jobs and empty jobs.
In mode 1, one Allgather exchanges those bounded summaries. Every rank resolves the small
rank-prefix summaries, then prepares its own job dictionaries independently.
This supports short outputs, marker-dependent tails, empty ranks and rotated
rank order; it does not require all tails to become literal during DISCOVER.
Mode 2 composes summaries in file order through an interval tree, then sends
each rank its actual incoming dictionary down the tree. There are `P-1`
compositions, `P-1` summary sends and `P-1` dictionary sends, followed by a final
dictionary broadcast. Empty ranks contribute identity summaries. Non-power-of-
two rank counts, short dictionaries and epoch-dependent rank rotation are
supported. The tree uses the decoder's private communicator, without new RMA
windows or a barrier at each level. Intermediate dictionaries are not gzip data
and do not change the decoded output.

CRC/ISIZE and member-boundary validation are unchanged. Summary storage and
scratch buffers are included in the decoder memory budget; compare initial job
counts and retries as well as time when testing different modes. The `history`
report decomposes the existing `window` timer, not an additional phase. MPI
time includes waiting for peers, not just network transfer. Decoder tests accept
`--history-mode 0|1|2` (`--parallel-history` remains an alias for 1).
`tests/test_mpi_gzip_history.cpp` checks tree incoming/final dictionaries against
a serial oracle, including empty jobs, short windows and communicator rotation.
`tests/benchmark_mpi_gzip_history.cpp` runs only the production decoder and its
complete gzip CRC/ISIZE checks, cycling the three modes without FASTQ parsing,
alignment or SAM output. The optional native test backend
`tests/backend_codec_single_cg.c` uses one 64-CPE group per MPI rank.
For a disk-independent benchmark, optionally link
`tests/benchmark_pread_cache.cpp` with `--wrap=pread --wrap=pread64` and set
`SWBWA_HISTORY_BENCH_CACHE=1`. An untimed warm-up captures this rank's compressed
reads (2 GiB by default, `SWBWA_HISTORY_CACHE_MIB=1..8192` for a test-only
override), then timed passes replay exactly those requests from
memory and reject cache misses. This does not change production input handling.
With a static MPI/C++ link, explicitly retaining `pthread_cond_broadcast`
(`-Wl,-u,pthread_cond_broadcast`) may be necessary for the standalone benchmark;
some toolchains otherwise leave libstdc++'s weak reference unresolved.

To benchmark the actual cross+pool dispatcher, use a fresh source/build
directory and keep the complete two-pass build:

```sh
MAKEFILES="$PWD/tests/benchmark_history.mk" \
    bash build.sh 8 EXEC_MODE=cgs_cross CPE_ALLOCATOR=pool USE_MPI=1 \
    HOST_MPE_THREADS=6 CPE_LDM_MODE=2 GZIP_INPUT=1 GZIP_OUTPUT=1 \
    LDFLAGS='tests/benchmark_history_runtime.o tests/benchmark_pread_cache.o -Wl,-u,pthread_cond_broadcast -Wl,--wrap=pread -Wl,--wrap=pread64'
```

This test-only make fragment replaces the main entry with the benchmark;
the resulting `SWBWA FILE.gz 3 384` is **not an alignment executable**.
It initializes the regular CPE pool and uses `backend_swbwa`, the cross
dispatcher, and the generated GOT/TLS relocation files. Submit from the build
directory with the normal cross-segment resource settings. Do not link a
standalone binary after layout discovery: that would invalidate the addresses.
Use a separate build directory for the regular application.

This is **bounded-window dynamic alignment scheduling**, not an arbitrary
gzip seek scheduler. Without prepare-ahead, window tails can still
wait for a slow reader; decoder metadata and dictionary exchanges remain
collective even when the optional output drain is disabled. Larger
windows trade memory for more scheduling freedom;
the chunk count cap is not an RSS cap. Besides existing decoder queues and
pipeline buffers, each rank retains at most `ready_chunks` source chunks (both
mates plus partitions), bounded individually by the formatter capacity. Packing
temporarily needs both the assembled chunks and window storage; originals are
released as they are copied, except the reserved chunk transferred to its owner.
Window capacity remains resident until close and retains its high-water size.
Packing, local copies and remote transfers can offset load-balance gains;
performance must be measured, not inferred from scheduling.

The `scheduler` report separates `refill_s`, `refill_wait_s`, `drain_wait_s`,
`claim_rma_s`, `fetch_s`, `pack_s`, and `local_copy_s`, with local/remote claims,
fetched bytes and peak packed ready-payload bytes (not peak RSS).
`direct_claims`/`direct_bytes` measure owner-local buffer transfers, not end-to-end
zero-copy decoding. `data_window_builds`, `window_capacity_bytes` and `window_s`
expose RMA allocation/reuse costs; `decode_drains` counts advances that actually
encountered pending output. `refill_s` includes packing, decoder drains and the
existing decode/framing/layout timers and must not be added to them.
`drain_wait_s` and `refill_wait_s` are waits, not Reader compute.
There are no extra per-chunk collectives in the claim path.

Host tests exercise the real decoder with a host codec launcher; they do not
replace native CPE/runtime validation:

```bash
python3 tests/test_mpi_gzip_decoder.py --work /tmp/swbwa-mpi-decoder-test
python3 tests/test_mpi_gzip_input.py --work /tmp/swbwa-mpi-input-test
python3 tests/test_mpi_gzip_input.py --input-mode dynamic --work /tmp/swbwa-mpi-dynamic-test
python3 tests/test_mpi_gzip_dynamic_pipeline.py --work /tmp/swbwa-mpi-dynamic-pipeline-test
python3 tests/test_mpi_gzip_pipeline_contract.py
```

Known native limitation: a 64-node gzip-input run with six `single + system`
MPI ranks per node (384 ranks total) has failed the `bns_fetch_seq` interval
assertion. Its cause remains unresolved; the host decoder and pipeline tests
do not reproduce or rule out this failure. Do not treat this configuration as
validated or the cross-segment placement changes as a fix for it.
