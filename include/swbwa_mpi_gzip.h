#ifndef SWBWA_MPI_GZIP_H
#define SWBWA_MPI_GZIP_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "swbwa_fastq_partition.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct swbwa_mpi_gzip_reader swbwa_mpi_gzip_reader;
typedef struct {
    char *data[2];
    int64_t size[2], chunk_id, first_record, records;
    swbwa_fastq_partition_t *partitions;
    unsigned partition_count;
} swbwa_mpi_gzip_batch;

/* Collective on WORLD, duplicated internally. MPI must already be initialized.
 * Uses the rapid decoder; no zlib fallback or preliminary whole-file scan.
 * path2 == NULL selects SE. batch_bytes is an R1 target, not a hard chunk size;
 * records/chunk is estimated once and then fixed. workers is also the CPE
 * formatter count (1..384).
 * Strict four LF-terminated lines/record, at most 1 MiB including LFs; no CR
 * or NUL. PE names are exact tokens before space/tab, without /1,/2 stripping.
 * The reader checks boundaries; the CPE formatter validates all records and
 * PE names using partitions. Fatal errors call MPI_Abort.
 * SWBWA_GZIP_DYNAMIC_TRANSFER selects staged (legacy redistribution), direct
 * (claim then gather fragments from decoder ranks), or local (direct with a
 * locality-matched first frontier). The latter two require PREPARE_AHEAD=0
 * and leave the decoder/history strategy unchanged. EARLY_FIRST=1 publishes
 * the first frontier before preparing its dynamically claimed suffix. They
 * are used only with dynamic input and more than one MPI rank.
 */
swbwa_mpi_gzip_reader *swbwa_mpi_gzip_open(const char *path1, const char *path2,
                                        size_t batch_bytes, unsigned workers, int debug);
/* STATIC: collective in the same order on every rank, READY may be empty.
 * DYNAMIC: independently claim increasing IDs from a bounded ready window;
 * refilling a window is collective, but the numbers of READY calls may differ.
 * Source-direct modes publish fragments instead of preassembled chunks. Bulk
 * FASTQ bytes move at most once across ranks, to the actual consumer; small
 * record-boundary carry transfers are separate. Delivered IDs remain increasing
 * on each rank, including the locality-matched initial frontier.
 * 1 = READY,
 * 0 = globally validated EOF; -1 is reserved for pre-MPI misuse only.
 * Each nonempty data[] is independently malloc-owned by the caller and has
 * a trailing NUL outside size[]. Records/first_record are zero-based logical
 * read ordinals (pair ordinals in PE), NOT compressed or decoded byte offsets.
 * Earlier batches remain provisional until successful EOF (gzip CRC checked).
 * partitions is separately malloc-owned. It must accompany data to the CPE
 * formatter and be freed after formatting; it has partition_count entries.
 */
int swbwa_mpi_gzip_next(swbwa_mpi_gzip_reader *, swbwa_mpi_gzip_batch *);
/* Call once per nonempty chunk, in local delivery order, after its output has
 * been consumed (Stage 3). Thread-safe; no MPI on the normal path. STATIC and
 * single-rank readers are no-ops. The optional DECODE_DRAIN counterfactual
 * waits for completions before new codec collectives. Does not free buffers. */
int swbwa_mpi_gzip_complete(swbwa_mpi_gzip_reader *, int64_t chunk_id);
/* Exact global count of nonempty chunks, available only after validated EOF. */
int64_t swbwa_mpi_gzip_chunk_count(const swbwa_mpi_gzip_reader *);
const char *swbwa_mpi_gzip_error(const swbwa_mpi_gzip_reader *);
/* Local report, no collective. Timings are rank-local seconds, bytes decoded.
 * Call from the next() consumer or after it has finished, never concurrently.
 * While background preparation is active only scheduler fields are reported;
 * full decode metrics become available after next() reaches EOF and joins it.
 * close is collective during normal operation and precedes MPI_Finalize.
 * It requires global EOF and all completion callbacks to have finished.
 * EOF alone need not drain in-flight output; premature close aborts.
 * The reader does not initialize/finalize the codec runtime or MPI.
 */
int swbwa_mpi_gzip_report(swbwa_mpi_gzip_reader *, FILE *);
void swbwa_mpi_gzip_close(swbwa_mpi_gzip_reader *);

#ifdef __cplusplus
}
#endif
#endif
