#ifndef SWBWA_GZIP_H
#define SWBWA_GZIP_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct swbwa_gzip_reader swbwa_gzip_reader;

typedef struct {
    uint64_t input_bytes, decoded_bytes, delivered_bytes;
    uint64_t jobs, accepted_jobs, fallback_bytes, capacity_retries;
    uint64_t source_read_bytes, prefetch_bytes, peak_symbol_bytes;
    uint64_t copy_bytes, copy_calls, refills;
    double read_s, kernel_s, queue_s, serial_s, line_index_s, copy_s;
    /* Rapid-engine phases. io_s measures foreground compressed-input reads
     * (including prefetch waits/copies), not isolated disk service time.
     * prefetch_io_s measures background pread; it overlaps foreground work.
     * zlib has no separate I/O decomposition: these fields remain zero and
     * its input I/O is included in serial_s. Do not sum inclusive phases. */
    double io_s, discover_run_s, resolve_run_s, repair_s;
    double input_prepare_s, task_build_s, prefetch_io_s;
    double discover_wait_s, resolve_wait_s;
} swbwa_gzip_input_metrics;

typedef struct {
    /* records counts SE reads or PE pairs, NOT the alignment n_processed. */
    uint64_t chunks, records, last_records;
    int paired;
    swbwa_gzip_input_metrics input[2];
} swbwa_gzip_metrics;

/* Gzip files only, four LF-terminated lines per FASTQ record. No full FASTQ
 * semantic/name/layout validation. Missing final LF and incomplete groups of
 * four lines are errors; bytes are never synthesized. Concatenated gzip
 * members and trailing zero padding are accepted; corruption is an error.
 *
 * path2 == NULL selects SE. chunk_bytes: 64 KiB..64 MiB; workers: 1..384.
 * backend must be "rapid" or "zlib". It is not read from the environment.
 * NULL on failure; error(NULL) gives this thread's latest open/API error.
 * Runtime init/finish and launch/alignment serialization belong to the caller.
 * Each reader is sequential, not thread-safe. No global read counters change. */
swbwa_gzip_reader *swbwa_gzip_open(const char *path1, const char *path2,
                                 size_t chunk_bytes, unsigned workers,
                                 const char *backend);

/* Set the decoded FASTQ batch target independently of compressed codec chunk.
 * Call before the first read_next. Default: 16 MiB of R1. Valid: 1..INT64_MAX-1
 * (also limited by size_t). Returns 0/-1; errors are sticky. Complete records
 * are returned up to this target; a single larger record may exceed it.
 * R1 selects the complete-record count; PE refills R2 until that same count
 * is ready. R2 byte lengths may be larger/smaller than the batch target. */
int swbwa_gzip_set_batch_bytes(swbwa_gzip_reader *reader, size_t batch_bytes);

/* 1: owned malloc buffers, 0: EOF, -1: sticky error. All four output pointers
 * are required and must be distinct. Outputs are reset to NULL/0 before work.
 * n1/n2 count bytes INCLUDING terminal LF, EXCLUDING the extra NUL sentinel.
 * Caller may mutate/free each buffer independently, even after reader close.
 * SE always sets *r2=NULL, *n2=0. PE emits equal record counts in input order,
 * not equal byte lengths; different total record counts are an error.
 * Drain to 0 to validate the entire input. Earlier chunks cannot be retracted
 * if later corruption, truncation or a mate-count mismatch is discovered. */
int swbwa_gzip_read_next(swbwa_gzip_reader *reader,
                        char **r1, int64_t *n1, char **r2, int64_t *n2);
int swbwa_gzip_next(swbwa_gzip_reader *reader,
                   char **r1, int64_t *n1, char **r2, int64_t *n2);
const char *swbwa_gzip_error(const swbwa_gzip_reader *reader);
void swbwa_gzip_close(swbwa_gzip_reader *reader);

/* Before handing the latest PE batch to the existing shared-offset CPE
 * formatter, call this with the actual CPE count (NOT the codec workers).
 * Buffers must still be unmodified. Checks equal total bytes and equal record
 * ranks at the formatter's R1-selected boundaries using newline metadata.
 * 0: compatible, -1: sticky error; O(cpe_count) boundary/page scans, no full
 * record walk. It does not require identical interior record byte lengths.
 * This check is mandatory for shared-offset formatting; read_next alone only
 * guarantees equal record counts. SE is a no-op. cpe_count: 1..384. */
int swbwa_gzip_check_pe_layout(swbwa_gzip_reader *reader,
                              const char *r1, int64_t n1,
                              const char *r2, int64_t n2,
                              unsigned cpe_count);

/* Snapshot/report: 0 success, -1 invalid argument or report I/O failure.
 * Metrics remain available after EOF/error until close. No implicit logging. */
int swbwa_gzip_get_metrics(const swbwa_gzip_reader *reader,
                          swbwa_gzip_metrics *metrics);
int swbwa_gzip_report(const swbwa_gzip_reader *reader, FILE *stream);

#ifdef __cplusplus
}
#endif
#endif
