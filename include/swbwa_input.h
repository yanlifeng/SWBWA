#ifndef SWBWA_INPUT_H
#define SWBWA_INPUT_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    unsigned int samples;
    double user_seconds, system_seconds;
    int64_t minor_faults, major_faults, voluntary_switches, involuntary_switches;
} swbwa_input_usage_t;

typedef struct {
    int fd;
    char *buffer;
    int64_t offset;
    size_t length;
} swbwa_input_slice_t;

typedef struct {
    uint64_t bytes, calls;
    double wall_seconds, syscall_seconds, max_seconds;
    int max_input;
    int64_t max_offset;
    uint64_t max_bytes;
    swbwa_input_usage_t usage;
} swbwa_input_result_t;

/* Parse SWBWA_INPUT_READERS: NULL selects all helpers, or fread for a one-MPE
 * build. Zero selects fread; invalid values return -1 with errno = EINVAL. */
int swbwa_input_reader_count(const char *value);

/* Thread-only CPU/resource accounting; samples == 0 means unavailable. */
void swbwa_input_usage_sample(swbwa_input_usage_t *usage);
void swbwa_input_usage_add(swbwa_input_usage_t *total,
                          const swbwa_input_usage_t *before,
                          const swbwa_input_usage_t *after);

/* Exact positioned reads into disjoint slices; helpers never allocate or call MPI.
 * The host worker pool must be initialized before requesting multiple readers. */
int swbwa_input_read(const swbwa_input_slice_t *inputs, int input_count,
                     int readers, int profile, swbwa_input_result_t *result);

#endif
