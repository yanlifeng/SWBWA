#ifndef SWBWA_MPI_DECODER_H
#define SWBWA_MPI_DECODER_H

#include <mpi.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct swbwa_mpi_decoder swbwa_mpi_decoder;

typedef struct {
    size_t chunk_bytes;          /* Compressed task stride; default 256 KiB. */
    size_t lookahead_bytes;      /* Bounded compressed block lookahead; 1 MiB. */
    size_t memory_bytes;         /* Adapter budget; 512 MiB per 64 workers. */
    size_t decoded_bytes;        /* Decoded budget; 128 MiB per 64 workers. */
    size_t initial_output_bytes; /* Initial decoded symbols per job; 2 MiB. */
    unsigned jobs_per_rank;      /* Maximum queued jobs, 1..384; matches workers. */
    unsigned workers;           /* Passed unchanged to codec launch; 1..384. */
} swbwa_mpi_decoder_config;

typedef struct {
    const unsigned char *data;
    size_t size;
    uint64_t global_offset;     /* Uncompressed stream offset, not a record ID. */
    uint64_t compressed_start_bit, compressed_end_bit;
    const uint32_t *line_counts; /* Newlines in each local 4096-byte page. */
    size_t line_count_pages;     /* Last page may be partial. No host rescan. */
    uint32_t crc32;
} swbwa_mpi_decoded_extent;

typedef struct {
    const unsigned char *bytes;  /* One contiguous rank-local decoded interval. */
    size_t size;
    uint64_t decoded_offset;     /* Offset of bytes[0]; ignore when size == 0. */
    uint64_t global_decoded_end; /* Committed global byte end after this epoch. */
    const swbwa_mpi_decoded_extent *extents;
    size_t count;                /* May be zero even when next returns 1. */
    uint64_t epoch;              /* Zero-based collective epoch. */
    uint64_t global_begin;       /* Global byte interval committed this call. */
    uint32_t stream_crc32;       /* CRC of [0, global_decoded_end). */
    int eof;                    /* Global, trailer-validated EOF. */
} swbwa_mpi_decoded_batch;

typedef struct {
    uint64_t epochs, discover_jobs, accepted_jobs, speculative_accepted;
    uint64_t discarded_jobs, capacity_retries, source_read_bytes, decoded_bytes;
    uint64_t peak_buffer_bytes;
    uint64_t initial_batch_jobs, max_discover_jobs, max_resolve_jobs;
    uint64_t known_tail_jobs, repaired_tail_jobs;
    uint64_t parallel_history_epochs, serial_history_epochs, tree_history_epochs;
    uint64_t history_mode;
    /* Accounted conservative scratch peak, including the H1 trace cache.
     * Tree payload/messages count point-to-point sends, excluding final Bcast;
     * summary bytes count locally constructed summaries, not MPI wire traffic. */
    uint64_t history_scratch_peak_bytes, history_summary_bytes;
    uint64_t history_tree_sent_bytes, history_tree_messages;
    uint64_t history_tree_compose_calls, history_tree_apply_calls;
    uint64_t prefetch_batches, prefetch_bytes, source_read_calls;
    double source_io_s, prefetch_io_s;
    uint64_t source_cache_bytes;
    double source_cache_s;
    int prefetch_cg;
    double io_s, discover_run_s, discover_wait_s, resolve_run_s, resolve_wait_s;
    /* Local job counts: discover = accepted + discarded + retry on success.
     * discarded: unselected jobs in committed attempts; retry: all jobs in
     * capacity-abandoned attempts. failed/full include both kinds of attempt.
     * discarded_output_bytes counts only successfully decoded, unselected
     * symbols in committed attempts, not partial output or retry attempts. */
    uint64_t retry_jobs, discover_failed_jobs, discover_output_full_jobs;
    uint64_t discarded_output_bytes;
    uint64_t gap_repair_jobs, gap_repair_ok;
    double gap_repair_s; /* Inclusive subset of next_s and backend/MPI timers. */
    /* Per-rank accumulated wall seconds, never summed across CPEs/ranks.
     * open_s includes open_io_s; io_s is foreground input time, including
     * waiting for prefetched bytes. source_io_s measures actual pread time;
     * prefetch_io_s is its overlapping background subset, not a next_s phase.
     * next_s includes retries and repeated EOF calls. Its disjoint phases
     * are io_s-open_io_s, prepare/cleanup, discover/resolve_wall, and the
     * metadata/window/commit/EOF timers. Backend run/wait are subsets of
     * discover/resolve_wall; retry_s is an inclusive subset of next_s.
     * MPI timers include time waiting for peers, not just network transfer. */
    double open_s, open_io_s, next_s, prepare_s, cleanup_s;
    double discover_wall_s, resolve_wall_s;
    /* Metadata packing, retry/chain planning vs Gather and plan broadcasts. */
    double metadata_local_s, metadata_mpi_s;
    /* Dictionary copies/tail marker repair/task selection vs broadcasts. */
    double window_local_s, window_mpi_s;
    /* Inclusive history wall; the remaining history timers are disjoint
     * subsets. MPI elapsed includes arrival/progress waits. Serial mode uses
     * window_mpi_s for its owner broadcasts; parallel_history counts mode 1. */
    double history_wall_s, history_summary_s, history_exchange_s;
    double history_linear_apply_s, history_job_dict_s;
    double history_tree_compose_s, history_tree_apply_s, history_tree_mpi_s;
    double history_final_sync_s;
    /* Extents/CRC/header validation (excluding pread) vs digest Gather,
     * commit Bcast and the final validation barrier, counted separately. */
    double commit_local_s, commit_mpi_s, commit_barrier_s, eof_barrier_s;
    double retry_s;
} swbwa_mpi_decoder_metrics;

/* Match the job and buffer budgets to the available CPEs. Allocations remain
 * demand-driven. Invalid worker counts are rejected by decoder_open. */
void swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *config, unsigned workers);

/* Collective on intracommunicator comm: duplicates it, never uses WORLD.
 * Every rank supplies identical options, path and gzip content. config=NULL
 * selects 384-worker defaults. Supports seekable, immutable gzip files, including
 * concatenated members. Every member's CRC/ISIZE is checked; trailing garbage
 * and padding are rejected. Inputs must remain unchanged during decoding;
 * filesystem timestamps and identity are not monitored. No pre-scan or
 * fallback. Empty members can yield a successful zero-byte epoch with eof=0;
 * all ranks must continue calling.
 * Budget covers adapter buffers, excluding MPI/runtime/codec worker storage.
 * open returns 0 on success, -1 on error and writes a bounded error message.
 * On ANY error the caller MUST promptly MPI_Abort its parent communicator:
 * peers may already be waiting in a different collective. Do not call close
 * or attempt recovery. No function here calls MPI_Abort, init, or finish.
 */
int swbwa_mpi_decoder_open(MPI_Comm comm, const char *path,
                          const swbwa_mpi_decoder_config *config,
                          swbwa_mpi_decoder **decoder, char *error, size_t error_size);

/* Collective, in identical order on all ranks. Returns 1 for a committed
 * epoch (including the last, possibly empty epoch with eof=1), 0 on subsequent
 * EOF calls, -1 on error. All ranks receive the same success/EOF status.
 * Each rank owns ONE contiguous interval, never interleaved with another
 * rank. extents are optional per-job newline indexes into the same bytes
 * allocation, in order; no FASTQ boundary is implied. Nonempty intervals
 * globally partition [global_begin, global_decoded_end). Rank order rotates.
 * Zero-byte DEFLATE jobs are committed but omitted from extents.
 * Returned read-only storage belongs to decoder and expires at the next
 * next/close. Copy before retaining across calls or formatting in place.
 * Earlier epochs remain provisional until a trailer-validated EOF is seen.
 * Invoke from a thread allowed by the application's MPI thread level.
 */
int swbwa_mpi_decoder_next(swbwa_mpi_decoder *decoder, swbwa_mpi_decoded_batch *batch);
const char *swbwa_mpi_decoder_error(const swbwa_mpi_decoder *decoder);
int swbwa_mpi_decoder_get_metrics(const swbwa_mpi_decoder *decoder,
                                swbwa_mpi_decoder_metrics *metrics);
/* Collective during normal operation only; call before MPI_Finalize.
 * Returns 0 and destroys decoder on success (NULL is also successful).
 * Returns -1 on failure, retaining decoder for error(); the caller MUST
 * immediately MPI_Abort, not retry close. Never ignores MPI_Comm_free errors.
 */
int swbwa_mpi_decoder_close(swbwa_mpi_decoder *decoder);

#ifdef __cplusplus
}
#endif
#endif
