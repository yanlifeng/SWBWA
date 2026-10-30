#ifndef SWBWA_OUTPUT_H
#define SWBWA_OUTPUT_H

#include <stddef.h>
#include <stdint.h>
#include "swbwa_config.h"
#include "bwa.h"

int swbwa_output_open(const char *path, int debug_enabled);
int swbwa_output_write(const void *data, size_t length);
/* One chunk owns one contiguous extent; single_ordered follows input IDs. */
int swbwa_output_write_chunk(int64_t id, int64_t start, int64_t end,
                             const bseq1_t *seqs, const int *sam_lengths,
                             int reads, int paired);
/* Producer-only, synchronous CPE compression; no output MPI or file operations.
 * Gzip preparation is enabled by default; SWBWA_GZIP_PREPARE_OUTPUT=0 disables it.
 * On success, *handle is NULL for disabled/plain output, otherwise an owned batch.
 * Returns 0 on success, -1/errno on failure with *handle left NULL.
 * At most SWBWA_PIPELINE_BUFFER_COUNT batches may be live; exhaustion is EAGAIN.
 * Reserve pipeline space before preparing. Do not hold the CPE runtime lock.
 * Encoder calls must stay on their runtime-owning thread, not arbitrary helpers. */
int swbwa_output_prepare_chunk(const bseq1_t *seqs, const int *sam_lengths,
                               int reads, void **handle);
/* Producer-only, after prepare and before queue publication. Optional ordered
 * gzip prefix-ahead mode enqueues only ID/size metadata, never SAM bytes.
 * Announced handles MUST be consumed by write_prepared_chunk in FIFO order;
 * they cannot be cancelled/released unwritten. Disabled modes are no-ops. */
int swbwa_output_publish_prepared_chunk(int64_t id, int reads, void *handle);
/* Writer-only. Consumes a non-NULL handle on success and after write failure;
 * an announced handle rejected before prefix consumption remains owned and the
 * pipeline must abort (a published extent cannot be cancelled). NULL falls back
 * to write_chunk. SAM is no longer needed for prepared real output, but RMA-only
 * validation still borrows seqs/lengths. Preserve FIFO chunk order. */
int swbwa_output_write_prepared_chunk(int64_t id, int64_t start, int64_t end,
                                      const bseq1_t *seqs, const int *sam_lengths,
                                      int reads, int paired, void *handle);
/* Return an unwritten handle (NULL is a no-op). No MPI; never use it again.
 * Drain/release all handles and join producer/writer before closing output;
 * close rejects outstanding prepared handles with EBUSY.
 * Legacy encoding APIs must not run concurrently with producer preparation. */
int swbwa_output_release_prepared_chunk(void *handle);
int swbwa_output_flush(void);
int swbwa_output_close(void);
const char *swbwa_output_name(void);

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
/* Online inputs publish monotonic chunk IDs without a preliminary count scan.
 * Finish is local, after validated input EOF and after the writer has drained. */
int swbwa_output_open_online(const char *path, int debug_enabled);
int swbwa_output_ordered_finish(int64_t chunks);
/* Reader-only: publish the zero-length prefix for an aligned empty chunk. */
int swbwa_output_ordered_skip(int64_t id);
#endif

#if SWBWA_OUTPUT_RMA_ONLY
int swbwa_output_begin_chunk(int64_t id, int64_t start, int64_t end,
                             int reads, int paired);
int swbwa_output_end_chunk(void);
#endif

#endif /* SWBWA_OUTPUT_H */
