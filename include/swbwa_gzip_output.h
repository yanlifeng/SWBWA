#ifndef SWBWA_GZIP_OUTPUT_H
#define SWBWA_GZIP_OUTPUT_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/uio.h>
#include "bwa.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct swbwa_gzip_output swbwa_gzip_output;
typedef struct swbwa_gzip_output_batch swbwa_gzip_output_batch;

swbwa_gzip_output *swbwa_gzip_output_create(unsigned workers);
/* One encoder caller; independently releasable buffers for a concurrent consumer.
 * EAGAIN means the bounded pool is full; EBUSY rejects concurrent encoders.
 * Input is borrowed only until return. No output-layer MPI/file operations;
 * the codec backend retains its existing runtime/progress contract. */
int swbwa_gzip_output_prepare(swbwa_gzip_output *writer, const bseq1_t *seqs,
                              const int *lengths, int reads,
                              swbwa_gzip_output_batch **batch);
int swbwa_gzip_output_batch_data(swbwa_gzip_output *writer,
                                 swbwa_gzip_output_batch *batch, int reads,
                                 const struct iovec **vectors, unsigned *count,
                                 uint64_t *raw_bytes, uint64_t *gzip_bytes);
int swbwa_gzip_output_release(swbwa_gzip_output *writer, swbwa_gzip_output_batch *batch);
/* Lifecycle check, after producer/consumer quiescence; excludes legacy borrowing. */
int swbwa_gzip_output_pending(swbwa_gzip_output *writer);
/* Returned iovecs borrow reusable compressed storage until the next encode.
 * SAM input is borrowed only during this synchronous call, never copied.
 * Do not overlap legacy borrowing with another encode/prepare call. */
int swbwa_gzip_output_encode(swbwa_gzip_output *writer, const bseq1_t *seqs,
                             const int *lengths, int reads,
                             const struct iovec **vectors, unsigned *count,
                             uint64_t *raw_bytes, uint64_t *gzip_bytes);
int swbwa_gzip_output_bytes(swbwa_gzip_output *writer, const void *data, size_t bytes,
                            const struct iovec **vectors, unsigned *count,
                            uint64_t *gzip_bytes);
void swbwa_gzip_output_report(const swbwa_gzip_output *writer, FILE *file,
                              int rank, int ranks);
void swbwa_gzip_output_destroy(swbwa_gzip_output *writer);

#ifdef __cplusplus
}
#endif

#endif
