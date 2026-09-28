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
int swbwa_output_flush(void);
int swbwa_output_close(void);
const char *swbwa_output_name(void);

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
/* Reader-only: publish the zero-length prefix for an aligned empty chunk. */
int swbwa_output_ordered_skip(int64_t id);
#endif

#if SWBWA_OUTPUT_RMA_ONLY
int swbwa_output_begin_chunk(int64_t id, int64_t start, int64_t end,
                             int reads, int paired);
int swbwa_output_end_chunk(void);
#endif

#endif /* SWBWA_OUTPUT_H */
