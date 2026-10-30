#ifndef SWBWA_CODEC_JOB_H
#define SWBWA_CODEC_JOB_H
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
enum { SWBWA_CODEC_DEFLATE = 1, SWBWA_CODEC_DISCOVER = 2,
       SWBWA_CODEC_RESOLVE = 3, SWBWA_CODEC_GZIP = 4 };
typedef struct {
    const unsigned char *data;
    size_t size;
} swbwa_codec_span;
typedef struct {
    const unsigned char *input;
    unsigned char *output;
    const unsigned char *dictionary;
    size_t input_size, output_capacity, output_size;
    uint32_t dictionary_size, crc;
    int operation, level, status;
    /* First-pass rapidgzip: bit offsets are relative to input. */
    uint16_t *symbols;
    size_t logical_size, start_bit, start_max_bit, stop_bit, end_bit;
    /* 0: search, 1: exact start with dictionary, 2: exact start with markers. */
    int known_start, final_block;
    unsigned char *tail;
    int tail_known;
    /* Distinguish output exhaustion from input bounds/invalid candidates. */
    int output_full;
    /* Optional newline counts for each 4096-byte decoded page. */
    uint32_t *line_counts;
    size_t state_bytes;
    int state_in_ldm;
    /* Borrowed SAM fragments: no packed input buffer is needed for compression. */
    const swbwa_codec_span *spans;
    unsigned span_count;
    size_t span_skip;
} swbwa_codec_job;
typedef struct {
    swbwa_codec_job *jobs;
    unsigned count, workers;
    double run_s, wait_s;
    unsigned char *workspaces;
    size_t workspace_bytes;
} swbwa_codec_batch;
void swbwa_codec_execute(swbwa_codec_job *job);
void swbwa_codec_gzip_worker(swbwa_codec_batch *batch, unsigned worker);
void swbwa_codec_launch(swbwa_codec_batch *batch);
void swbwa_codec_init(void);
void swbwa_codec_finish(void);
#ifdef __cplusplus
}
#endif
#endif
