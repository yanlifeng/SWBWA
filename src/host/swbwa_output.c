#include "swbwa_config.h"
#include "swbwa_mpi.h"
#include "swbwa_output.h"
#include "swbwa_discard_digest.h"
#include "swbwa_host_workers.h"
#include "swbwa_sam_md5.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

#if SWBWA_OUTPUT_RMA_ONLY
enum { OUTPUT_SAMPLE_READS = 100 };
typedef struct {
    int64_t id, start, end;
    int expected_reads, paired;
    uint64_t reads, bytes, sample_reads, sample_bytes, sum, xor;
} swbwa_output_sample_t;
#endif

typedef struct {
    int64_t id, start, end;
    uint64_t offset, bytes;
    int reads;
} swbwa_output_extent_t;

typedef struct {
    unsigned char *buffer;
    size_t used;
    size_t capacity;
    uint64_t write_calls;
    uint64_t submitted_bytes;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    uint64_t hash_sum;
    uint64_t hash_xor;
#endif
    uint64_t buffered_flush_calls;
    uint64_t buffered_flush_bytes;
    uint64_t direct_write_calls;
    uint64_t direct_write_bytes;
    uint64_t posix_write_calls;
    uint64_t posix_write_bytes;
    uint64_t reservation_calls;
    uint64_t reserved_bytes;
    uint64_t pwrite_calls;
    uint64_t pwrite_bytes;
    uint64_t last_offset;
    size_t *record_lengths;
    size_t record_capacity;
    swbwa_output_extent_t *extents;
    size_t extent_count, extent_capacity;
    uint64_t chunk_calls;
    double chunk_measure_seconds, chunk_pack_seconds;
    double buffered_flush_seconds;
    double posix_write_seconds;
    double fetch_and_op_seconds;
    double win_flush_seconds;
    double reservation_total_seconds;
    double pwrite_seconds;
    double win_unlock_all_seconds;
    double win_free_seconds;
    double file_close_seconds;
    int debug_enabled;
    int opened;
    int fd;
    int owns_fd;
    int stream_md5;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    int discard_hash_enabled;
#endif
    char *name;
#if SWBWA_OUTPUT_RMA_ONLY
    int chunk_active;
    swbwa_output_sample_t chunk;
    swbwa_output_sample_t *samples;
    size_t sample_count, sample_capacity;
#endif
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    MPI_Win offset_window;
    uint64_t *offset_base;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    int64_t ordered_count, last_chunk;
    uint64_t ordered_empty_chunks, ordered_polls;
    double ordered_wait_seconds;
#endif
#endif
} swbwa_output_state_t;

static swbwa_output_state_t output_state;

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
static int discard_hash_enabled(void)
{
    const char *value = getenv("SWBWA_DISCARD_HASH");

    if (value == NULL || *value == '\0' || strcmp(value, "1") == 0)
        return 1;
    if (strcmp(value, "0") == 0) return 0;

    if (swbwa_mpi_is_root())
        fprintf(stderr,
                "[E::output] SWBWA_DISCARD_HASH must be 0 or 1, got '%s'\n",
                value);
    errno = EINVAL;
    return -1;
}
#endif

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD || SWBWA_OUTPUT_RMA_ONLY
static inline uint64_t hash_sam_word(const unsigned char *cursor)
{
    const uint64_t multiplier = UINT64_C(0xc6a4a7935bd1e995);
    uint64_t word;

    memcpy(&word, cursor, sizeof(word));
    word *= multiplier;
    word ^= word >> 47;
    return word * multiplier;
}

static uint64_t hash_sam_record(const void *data, size_t length)
{
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
    return swbwa_digest_hash_blob(data, length);
#else
    static const uint64_t multiplier = UINT64_C(0xc6a4a7935bd1e995);
    static const unsigned int shift = 47;
    const unsigned char *cursor = data;
    size_t remaining = length;
    uint64_t hash = UINT64_C(0x9e3779b97f4a7c15) ^
                    ((uint64_t)length * multiplier);

    /* Independent word mixing can overlap; the record recurrence is unchanged. */
    while (remaining >= 4 * sizeof(uint64_t)) {
        uint64_t word0 = hash_sam_word(cursor);
        uint64_t word1 = hash_sam_word(cursor + 8);
        uint64_t word2 = hash_sam_word(cursor + 16);
        uint64_t word3 = hash_sam_word(cursor + 24);

        hash = (hash ^ word0) * multiplier;
        hash = (hash ^ word1) * multiplier;
        hash = (hash ^ word2) * multiplier;
        hash = (hash ^ word3) * multiplier;
        cursor += 4 * sizeof(uint64_t);
        remaining -= 4 * sizeof(uint64_t);
    }
    while (remaining >= sizeof(uint64_t)) {
        hash = (hash ^ hash_sam_word(cursor)) * multiplier;
        cursor += sizeof(uint64_t);
        remaining -= sizeof(uint64_t);
    }

    switch (remaining) {
    case 7: hash ^= (uint64_t)cursor[6] << 48; /* fall through */
    case 6: hash ^= (uint64_t)cursor[5] << 40; /* fall through */
    case 5: hash ^= (uint64_t)cursor[4] << 32; /* fall through */
    case 4: hash ^= (uint64_t)cursor[3] << 24; /* fall through */
    case 3: hash ^= (uint64_t)cursor[2] << 16; /* fall through */
    case 2: hash ^= (uint64_t)cursor[1] << 8;  /* fall through */
    case 1:
        hash ^= (uint64_t)cursor[0];
        hash *= multiplier;
        break;
    default:
        break;
    }

    hash ^= hash >> shift;
    hash *= multiplier;
    hash ^= hash >> shift;
    return hash;
#endif
}
#endif

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
static void print_discard_hash(void)
{
    fprintf(stderr,
            "[SWBWA output hash rank %06d/%06d] calls=%" PRIu64
            " bytes=%" PRIu64 " sum=0x%016" PRIx64
            " xor=0x%016" PRIx64 " enabled=%d hash_prefix_bytes=%llu\n",
            swbwa_mpi_rank(), swbwa_mpi_size(),
            output_state.write_calls, output_state.submitted_bytes,
            output_state.hash_sum, output_state.hash_xor,
            output_state.discard_hash_enabled,
            (unsigned long long)SWBWA_DISCARD_HASH_BYTES);
    fflush(stderr);
}
#endif

#if SWBWA_OUTPUT_RMA_ONLY
int swbwa_output_begin_chunk(int64_t id, int64_t start, int64_t end,
                             int reads, int paired)
{
    if (!output_state.opened || output_state.chunk_active || id < 0 ||
        start < 0 || end < start || reads < 0 ||
        (paired != 0 && paired != 1) || (paired && (reads & 1))) {
        errno = EINVAL;
        return -1;
    }
    memset(&output_state.chunk, 0, sizeof(output_state.chunk));
    output_state.chunk.id = id;
    output_state.chunk.start = start;
    output_state.chunk.end = end;
    output_state.chunk.expected_reads = reads;
    output_state.chunk.paired = paired;
    output_state.chunk_active = 1;
    return 0;
}

int swbwa_output_end_chunk(void)
{
    if (!output_state.opened || !output_state.chunk_active ||
        output_state.chunk.reads != (uint64_t)output_state.chunk.expected_reads) {
        errno = EINVAL;
        return -1;
    }
    if (output_state.sample_count == output_state.sample_capacity) {
        size_t capacity = output_state.sample_capacity ?
                          output_state.sample_capacity * 2 : 16;
        swbwa_output_sample_t *samples;

        if (capacity < output_state.sample_capacity ||
            capacity > SIZE_MAX / sizeof(*samples)) {
            errno = EOVERFLOW;
            return -1;
        }
        samples = realloc(output_state.samples, capacity * sizeof(*samples));
        if (samples == NULL) return -1;
        output_state.samples = samples;
        output_state.sample_capacity = capacity;
    }
    output_state.samples[output_state.sample_count++] = output_state.chunk;
    output_state.chunk_active = 0;
    /* Chunk boundaries must not change the normal output flush schedule. */
    return 0;
}

static void print_output_samples(void)
{
    size_t i;

    for (i = 0; i < output_state.sample_count; ++i) {
        const swbwa_output_sample_t *s = &output_state.samples[i];
        fprintf(stderr,
                "[SWBWA chunk sample rank %06d/%06d] chunk=%" PRId64
                " start=%" PRId64 " end=%" PRId64 " paired=%d"
                " reads=%" PRIu64 " sam_bytes=%" PRIu64
                " sample_reads=%" PRIu64 " sample_bytes=%" PRIu64
                " sum=0x%016" PRIx64 " xor=0x%016" PRIx64 "\n",
                swbwa_mpi_rank(), swbwa_mpi_size(), s->id, s->start, s->end,
                s->paired, s->reads, s->bytes, s->sample_reads, s->sample_bytes,
                s->sum, s->xor);
    }
    fprintf(stderr,
            "[SWBWA output RMA-only rank %06d/%06d] chunks=%zu"
            " reads=%" PRIu64 " sam_bytes=%" PRIu64
            " reserved_bytes=%" PRIu64 " reservations=%" PRIu64
            " disk_bytes=%" PRIu64 " sample_limit=%d\n",
            swbwa_mpi_rank(), swbwa_mpi_size(), output_state.sample_count,
            output_state.write_calls, output_state.submitted_bytes,
            output_state.reserved_bytes, output_state.reservation_calls,
            output_state.pwrite_bytes, OUTPUT_SAMPLE_READS);
}
#endif

static double output_debug_now(void)
{
    struct timespec now;

    if (!output_state.debug_enabled) return 0.0;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec * 1.0e-9;
}

static char *make_split_name(const char *path, int rank)
{
    static const char suffix[] = ".sam";
    size_t length = strlen(path);
    size_t stem_length = length;
    size_t capacity;
    char *name;

    if (length >= sizeof(suffix) - 1 &&
        strcmp(path + length - (sizeof(suffix) - 1), suffix) == 0)
        stem_length -= sizeof(suffix) - 1;
    capacity = stem_length + sizeof(".rank000000.sam") + 16;
    name = malloc(capacity);
    if (name == NULL) return NULL;
    snprintf(name, capacity, "%.*s.rank%06d.sam", (int)stem_length, path, rank);
    return name;
}

static int write_all(int fd, const unsigned char *data, size_t length)
{
    if (output_state.stream_md5)
        return swbwa_sam_md5_update(data, length);
    while (length > 0) {
        double start = output_debug_now();
        ssize_t written = write(fd, data, length);

        if (output_state.debug_enabled) {
            ++output_state.posix_write_calls;
            output_state.posix_write_seconds += output_debug_now() - start;
        }
        if (written < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (written == 0) {
            errno = EIO;
            return -1;
        }
        if (output_state.debug_enabled)
            output_state.posix_write_bytes += (uint64_t)written;
        data += written;
        length -= (size_t)written;
    }
    return 0;
}

#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
static int pwrite_all(int fd, const unsigned char *data, size_t length,
                      uint64_t offset)
{
    while (length > 0) {
        size_t chunk = length > (size_t)INT_MAX ? (size_t)INT_MAX : length;
        double start = output_debug_now();
        ssize_t written = pwrite(fd, data, chunk, (off_t)offset);

        if (output_state.debug_enabled) {
            ++output_state.pwrite_calls;
            output_state.pwrite_seconds += output_debug_now() - start;
        }
        if (written < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (written == 0) {
            errno = EIO;
            return -1;
        }
        if (output_state.debug_enabled)
            output_state.pwrite_bytes += (uint64_t)written;
        data += written;
        length -= (size_t)written;
        offset += (uint64_t)written;
    }
    return 0;
}

static int mpi_check(int result, const char *operation)
{
    char error[MPI_MAX_ERROR_STRING];
    int length = 0;

    if (result == MPI_SUCCESS) return 0;
    errno = EIO;
    MPI_Error_string(result, error, &length);
    fprintf(stderr, "[E::MPI output rank %d] %s failed: %.*s\n",
            swbwa_mpi_rank(), operation, length, error);
    return -1;
}

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
/* Global input tickets make each rank's FIFO monotonic. Therefore the first
 * unpublished prefix cannot be blocked behind a later chunk on that rank. */
static int reserve_ordered(int64_t id, uint64_t bytes, uint64_t *offset,
                           int record_stats)
{
    uint64_t next, polls = 0;
    double start = output_debug_now();
    int result;

    if (id < 0 || id >= output_state.ordered_count || bytes > INT64_MAX) {
        errno = EINVAL;
        return -1;
    }
    do {
        result = MPI_Fetch_and_op(NULL, offset, MPI_UINT64_T, 0, (MPI_Aint)id,
                                  MPI_NO_OP, output_state.offset_window);
        if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
        if (mpi_check(result, "read ordered prefix") != 0) return -1;
        ++polls;
        if (*offset == UINT64_MAX) sched_yield();
    } while (*offset == UINT64_MAX);
    if (*offset > (uint64_t)INT64_MAX - bytes) {
        errno = EOVERFLOW;
        return -1;
    }
    next = *offset + bytes;
    /* Accumulate/NO_OP provides atomic access to each prefix slot; ordinary
     * concurrent Put/Get would not provide this atomicity contract. */
    result = MPI_Accumulate(&next, 1, MPI_UINT64_T, 0, (MPI_Aint)(id + 1),
                            1, MPI_UINT64_T, MPI_REPLACE, output_state.offset_window);
    if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
    if (mpi_check(result, "publish ordered prefix") != 0) return -1;
    if (record_stats) {
        double seconds = output_debug_now() - start;
        ++output_state.reservation_calls;
        output_state.reserved_bytes += bytes;
        output_state.ordered_polls += polls;
        output_state.ordered_wait_seconds += seconds;
        output_state.reservation_total_seconds += seconds;
        output_state.last_offset = *offset;
    }
    return 0;
}

int swbwa_output_ordered_skip(int64_t id)
{
    uint64_t offset;
    if (!output_state.opened) {
        errno = EINVAL;
        return -1;
    }
    /* Only the reader updates this counter; do not touch writer statistics. */
    if (reserve_ordered(id, 0, &offset, 0) != 0) return -1;
    ++output_state.ordered_empty_chunks;
    return 0;
}
#endif

static int write_single_unordered(const unsigned char *data, size_t length)
{
    uint64_t increment;
    uint64_t offset;
    int result;
    double reservation_start;
    double operation_start;

    if (length == 0) return 0;
    increment = (uint64_t)length;
    if (increment > (uint64_t)INT64_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    reservation_start = output_debug_now();
    if (output_state.debug_enabled) ++output_state.reservation_calls;
    operation_start = output_debug_now();

    result = MPI_Fetch_and_op(&increment, &offset, MPI_UINT64_T, 0, 0,
                              MPI_SUM, output_state.offset_window);
    if (output_state.debug_enabled) {
        double operation_end = output_debug_now();

        output_state.fetch_and_op_seconds += operation_end - operation_start;
        operation_start = operation_end;
    }
    if (result == MPI_SUCCESS) {
        result = MPI_Win_flush(0, output_state.offset_window);
        if (output_state.debug_enabled)
            output_state.win_flush_seconds +=
                output_debug_now() - operation_start;
    }
    if (output_state.debug_enabled)
        output_state.reservation_total_seconds +=
            output_debug_now() - reservation_start;

    if (mpi_check(result, "MPI output offset reservation") != 0)
        return -1;
    if (offset > (uint64_t)INT64_MAX - increment) {
        errno = EOVERFLOW;
        return -1;
    }
    output_state.reserved_bytes += increment;
    output_state.last_offset = offset;

    /*
     * RMA gives every rank a disjoint file extent.  POSIX pwrite keeps those
     * writes parallel without routing bulk I/O through Sunway MPI progress.
     */
#if SWBWA_OUTPUT_RMA_ONLY
    (void)data;
    return 0;
#else
    return pwrite_all(output_state.fd, data, length, offset);
#endif
}

static int flush_single_unordered(void)
{
    if (write_single_unordered(output_state.buffer, output_state.used) != 0)
        return -1;
    output_state.used = 0;
    return 0;
}
#endif

int swbwa_output_open(const char *path, int debug_enabled)
{
    const char *md5_option = getenv("SWBWA_OUTPUT_MD5");
    int stream_md5 = md5_option != NULL && strcmp(md5_option, "1") == 0;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    size_t capacity = 0;
#else
    size_t capacity = (size_t)SWBWA_OUTPUT_BUFFER_BYTES;
#endif

    if (output_state.opened) {
        errno = EALREADY;
        return -1;
    }
    if (md5_option && strcmp(md5_option, "0") && strcmp(md5_option, "1")) {
        fprintf(stderr, "[E::output] SWBWA_OUTPUT_MD5 must be 0 or 1\n");
        errno = EINVAL;
        return -1;
    }
#if SWBWA_USE_MPI || SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD || SWBWA_OUTPUT_RMA_ONLY
    if (stream_md5) {
        fprintf(stderr, "[E::output] streaming MD5 requires non-MPI file output\n");
        errno = ENOTSUP;
        return -1;
    }
#endif
#if !(SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)
    if (capacity == 0 || capacity > INT_MAX) {
        errno = EINVAL;
        return -1;
    }
#endif
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    output_state.capacity = capacity;
    output_state.debug_enabled = debug_enabled != 0 || SWBWA_OUTPUT_RMA_ONLY;
#if !(SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)
    output_state.buffer = malloc(capacity);
    if (output_state.buffer == NULL) return -1;
#endif

    if (stream_md5) {
        output_state.name = strdup("(ordered SAM streaming MD5; no file)");
        if (!output_state.name || swbwa_sam_md5_open() != 0) goto fail;
        output_state.stream_md5 = 1;
        output_state.opened = 1;
        return 0;
    }

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    (void)path;
    output_state.discard_hash_enabled = discard_hash_enabled();
    if (output_state.discard_hash_enabled < 0) goto fail;
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
    {
        const uint16_t endian_probe = 1;
        if (*(const unsigned char *)&endian_probe != 1) {
            errno = EINVAL;
            goto fail;
        }
    }
    fprintf(stderr, "[SWBWA discard digest] hash_scope=generated_sam final_sam_copy=0 metadata_bytes=%d\n",
            SWBWA_DIGEST_RECORD_BYTES);
#endif
    output_state.name = strdup("(discard)");
    if (output_state.name == NULL) goto fail;
#elif SWBWA_USE_MPI
    if (path == NULL || strcmp(path, "-") == 0) {
        errno = EINVAL;
        goto fail;
    }
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SPLIT
    output_state.name = make_split_name(path, swbwa_mpi_rank());
    if (output_state.name == NULL) goto fail;
    output_state.fd = open(output_state.name, O_CREAT | O_WRONLY | O_TRUNC, 0666);
    if (output_state.fd < 0) goto fail;
    output_state.owns_fd = 1;
#elif SWBWA_OUTPUT_SINGLE_FILE
#if !SWBWA_OUTPUT_RMA_ONLY
    int local_open_failed;
    int local_open_errno;
    int any_open_failed;
    int open_flags = O_CREAT | O_WRONLY;
#endif

    if (sizeof(off_t) < sizeof(int64_t)) {
        errno = EOVERFLOW;
        goto fail;
    }
    output_state.name = strdup(path);
    if (output_state.name == NULL) goto fail;
#if !SWBWA_OUTPUT_RMA_ONLY
    if (swbwa_mpi_rank() == 0) open_flags |= O_TRUNC;
    output_state.fd = open(output_state.name, open_flags, 0666);
    output_state.owns_fd = output_state.fd >= 0;
    local_open_failed = output_state.fd < 0;
    local_open_errno = local_open_failed ? errno : 0;
    if (mpi_check(MPI_Allreduce(&local_open_failed, &any_open_failed, 1,
                                MPI_INT, MPI_MAX, MPI_COMM_WORLD),
                  "MPI_Allreduce after output open") != 0)
        goto fail;
    if (any_open_failed) {
        if (local_open_failed)
            fprintf(stderr, "[E::output rank %d] failed to open %s: %s\n",
                    swbwa_mpi_rank(), output_state.name,
                    strerror(local_open_errno));
        errno = EIO;
        goto fail;
    }
#endif
    size_t window_bytes = sizeof(uint64_t);
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    output_state.ordered_count = swbwa_mpi_fastq_scheduler_chunk_count();
    output_state.last_chunk = -1;
    if (output_state.ordered_count < 0 ||
        (uint64_t)output_state.ordered_count >= (uint64_t)PTRDIFF_MAX / sizeof(uint64_t) ||
        strcmp(swbwa_mpi_fastq_scheduler_ticket_mode(), "global") != 0) {
        errno = EINVAL;
        goto fail;
    }
    window_bytes = ((size_t)output_state.ordered_count + 1) * sizeof(uint64_t);
#endif
    if (mpi_check(MPI_Win_allocate(swbwa_mpi_rank() == 0 ? window_bytes : 0,
                                   sizeof(uint64_t), MPI_INFO_NULL, MPI_COMM_WORLD,
                                   &output_state.offset_base,
                                   &output_state.offset_window),
                  "MPI_Win_allocate") != 0)
        goto fail;
    if (mpi_check(MPI_Win_lock_all(0, output_state.offset_window),
                  "MPI_Win_lock_all") != 0)
        goto fail_window;
    if (swbwa_mpi_rank() == 0) {
        memset(output_state.offset_base, 0xff, window_bytes);
        *output_state.offset_base = 0;
    }
    if (mpi_check(MPI_Win_sync(output_state.offset_window), "MPI_Win_sync") != 0)
        goto fail_locked_window;
    if (mpi_check(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier") != 0)
        goto fail_locked_window;
#endif
#else
    if (path == NULL || strcmp(path, "-") == 0) {
        output_state.fd = STDOUT_FILENO;
        output_state.name = strdup("stdout");
    } else {
        output_state.name = strdup(path);
        if (output_state.name == NULL) goto fail;
        output_state.fd = open(path, O_CREAT | O_WRONLY | O_TRUNC, 0666);
        if (output_state.fd < 0) goto fail;
        output_state.owns_fd = 1;
    }
    if (output_state.name == NULL) goto fail;
#endif

    output_state.opened = 1;
    return 0;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
fail_locked_window:
    MPI_Win_unlock_all(output_state.offset_window);
fail_window:
    MPI_Win_free(&output_state.offset_window);
#endif
fail:
    if (output_state.owns_fd && output_state.fd >= 0) close(output_state.fd);
    free(output_state.name);
    free(output_state.buffer);
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    return -1;
}

#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
int swbwa_output_discard_hash_active(void)
{
    if (!output_state.opened) {
        errno = EINVAL;
        return -1;
    }
    return output_state.discard_hash_enabled;
}

int swbwa_output_write_digest(const void *record)
{
    const unsigned char *p = (const unsigned char *)record;
    uint64_t length, hash, state, expected;

    if (!output_state.opened || p == NULL) {
        errno = EINVAL;
        return -1;
    }
    length = swbwa_digest_load64(p);
    hash = swbwa_digest_load64(p + 8);
    state = swbwa_digest_load64(p + 16);
    expected = output_state.discard_hash_enabled ? SWBWA_DIGEST_HASH_READY : SWBWA_DIGEST_COUNT_READY;
    if (state != expected) {
        errno = EINVAL;
        return -1;
    }
    if (length == 0) return 0;
    ++output_state.write_calls;
    output_state.submitted_bytes += length;
    if (output_state.discard_hash_enabled) {
        output_state.hash_sum += hash;
        output_state.hash_xor ^= hash;
    }
    return 0;
}
#endif

int swbwa_output_write(const void *data, size_t length)
{
    const unsigned char *source = data;

    if (!output_state.opened || (data == NULL && length != 0)) {
        errno = EINVAL;
        return -1;
    }
    if (length == 0) return 0;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    /* A streaming write has no input chunk identity and cannot be ordered. */
    errno = EINVAL;
    return -1;
#endif
#if SWBWA_OUTPUT_RMA_ONLY
    {
        swbwa_output_sample_t *s = &output_state.chunk;
        if (!output_state.chunk_active || s->reads >= (uint64_t)s->expected_reads) {
            errno = EINVAL;
            return -1;
        }
        if (s->reads < OUTPUT_SAMPLE_READS) {
            uint64_t hash = hash_sam_record(data, length);
            ++s->sample_reads;
            s->sample_bytes += (uint64_t)length;
            s->sum += hash;
            s->xor ^= hash;
        }
        ++s->reads;
        s->bytes += (uint64_t)length;
    }
#endif
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    {
        ++output_state.write_calls;
        output_state.submitted_bytes += (uint64_t)length;
        if (output_state.discard_hash_enabled) {
            size_t hash_length = length;
            uint64_t hash;

#if SWBWA_DISCARD_HASH_BYTES > 0
            if (hash_length > SWBWA_DISCARD_HASH_BYTES)
                hash_length = SWBWA_DISCARD_HASH_BYTES;
#endif
            hash = hash_sam_record(source, hash_length);
            output_state.hash_sum += hash;
            output_state.hash_xor ^= hash;
        }
    }
    return 0;
#else
    if (output_state.debug_enabled) {
        ++output_state.write_calls;
        output_state.submitted_bytes += (uint64_t)length;
    }
#endif
    if (length > output_state.capacity - output_state.used &&
        swbwa_output_flush() != 0)
        return -1;
    if (length > output_state.capacity) {
        if (output_state.debug_enabled) {
            ++output_state.direct_write_calls;
            output_state.direct_write_bytes += (uint64_t)length;
        }
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
        return write_single_unordered(source, length);
#else
        return write_all(output_state.fd, source, length);
#endif
    }
    memcpy(output_state.buffer + output_state.used, source, length);
    output_state.used += length;
    if (output_state.used == output_state.capacity)
        return swbwa_output_flush();
    return 0;
}

#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
typedef struct {
    const bseq1_t *seqs;
    const int *known_lengths;
    size_t *lengths;
    unsigned char *packed;
    size_t sums[SWBWA_HOST_MPE_THREADS];
    size_t offsets[SWBWA_HOST_MPE_THREADS];
    int errors[SWBWA_HOST_MPE_THREADS];
    int reads;
} output_pack_t;

static void measure_chunk(void *opaque, int worker, int workers)
{
    output_pack_t *task = opaque;
    size_t begin = (size_t)task->reads * worker / workers;
    size_t end = (size_t)task->reads * (worker + 1) / workers;
    size_t total = 0, i;

    for (i = begin; i < end; ++i) {
        size_t length;
        if (task->seqs[i].sam == NULL) {
            task->errors[worker] = EINVAL;
            return;
        }
        if (task->known_lengths != NULL && task->known_lengths[i] < 0) {
            task->errors[worker] = EINVAL;
            return;
        }
        length = task->known_lengths != NULL ? (size_t)task->known_lengths[i]
                                             : strlen(task->seqs[i].sam);
        if (length > SIZE_MAX - total) {
            task->errors[worker] = EOVERFLOW;
            return;
        }
        task->lengths[i] = length;
        total += length;
    }
    task->sums[worker] = total;
}

static void pack_chunk(void *opaque, int worker, int workers)
{
    output_pack_t *task = opaque;
    size_t begin = (size_t)task->reads * worker / workers;
    size_t end = (size_t)task->reads * (worker + 1) / workers;
    size_t pos = task->offsets[worker], i;

    for (i = begin; i < end; ++i) {
        memcpy(task->packed + pos, task->seqs[i].sam, task->lengths[i]);
        pos += task->lengths[i];
    }
}
#endif

int swbwa_output_write_chunk(int64_t id, int64_t start, int64_t end,
                             const bseq1_t *seqs, const int *sam_lengths,
                             int reads, int paired)
{
    int i;
    if (!output_state.opened || id < 0 || start < 0 || end < start ||
        reads < 0 || (reads && seqs == NULL) ||
        (paired != 0 && paired != 1) || (paired && (reads & 1))) {
        errno = EINVAL;
        return -1;
    }
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    {
        output_pack_t task;
        size_t bytes = 0;
        double clock_start;

#if SWBWA_OUTPUT_RMA_ONLY
        if (output_state.chunk_active) {
            errno = EINVAL;
            return -1;
        }
#endif
        /* Drain legacy streaming output before giving this chunk its own extent. */
        if (swbwa_output_flush() != 0) return -1;
        if ((size_t)reads > output_state.record_capacity) {
            size_t *lengths;
            if ((size_t)reads > SIZE_MAX / sizeof(*lengths)) {
                errno = EOVERFLOW;
                return -1;
            }
            lengths = realloc(output_state.record_lengths, (size_t)reads * sizeof(*lengths));
            if (lengths == NULL) return -1;
            output_state.record_lengths = lengths;
            output_state.record_capacity = (size_t)reads;
        }
        memset(&task, 0, sizeof(task));
        task.seqs = seqs;
        task.known_lengths = sam_lengths;
        task.lengths = output_state.record_lengths;
        task.reads = reads;
        clock_start = output_debug_now();
        swbwa_host_workers_run_ready(measure_chunk, &task);
        for (i = 0; i < SWBWA_HOST_MPE_THREADS; ++i) {
            if (task.errors[i] || task.sums[i] > SIZE_MAX - bytes) {
                errno = task.errors[i] ? task.errors[i] : EOVERFLOW;
                return -1;
            }
            task.offsets[i] = bytes;
            bytes += task.sums[i];
        }
        output_state.chunk_measure_seconds += output_debug_now() - clock_start;
        if ((uint64_t)bytes > INT64_MAX) {
            errno = EOVERFLOW;
            return -1;
        }
        if (bytes > output_state.capacity) {
            unsigned char *buffer = realloc(output_state.buffer, bytes);
            if (buffer == NULL) return -1;
            output_state.buffer = buffer;
            output_state.capacity = bytes;
        }
        if (output_state.debug_enabled && output_state.extent_count == output_state.extent_capacity) {
            size_t capacity = output_state.extent_capacity ? output_state.extent_capacity * 2 : 16;
            swbwa_output_extent_t *extents;
            if (capacity < output_state.extent_capacity || capacity > SIZE_MAX / sizeof(*extents)) {
                errno = EOVERFLOW;
                return -1;
            }
            extents = realloc(output_state.extents, capacity * sizeof(*extents));
            if (extents == NULL) return -1;
            output_state.extents = extents;
            output_state.extent_capacity = capacity;
        }
        task.packed = output_state.buffer;
        clock_start = output_debug_now();
        swbwa_host_workers_run_ready(pack_chunk, &task);
        output_state.chunk_pack_seconds += output_debug_now() - clock_start;
#if SWBWA_OUTPUT_RMA_ONLY
        if (swbwa_output_begin_chunk(id, start, end, reads, paired) != 0) return -1;
        for (i = 0; i < reads && i < OUTPUT_SAMPLE_READS; ++i) {
            uint64_t hash = hash_sam_record(seqs[i].sam, task.lengths[i]);
            ++output_state.chunk.sample_reads;
            output_state.chunk.sample_bytes += task.lengths[i];
            output_state.chunk.sum += hash;
            output_state.chunk.xor ^= hash;
        }
        output_state.chunk.reads = (uint64_t)reads;
        output_state.chunk.bytes = bytes;
#endif
        ++output_state.chunk_calls;
        output_state.write_calls += (uint64_t)reads;
        output_state.submitted_bytes += bytes;
        output_state.used = bytes;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        uint64_t offset;
        if (id <= output_state.last_chunk) {
            errno = EINVAL;
            return -1;
        }
        if (reserve_ordered(id, bytes, &offset, 1) != 0) return -1;
#if !SWBWA_OUTPUT_RMA_ONLY
        if (pwrite_all(output_state.fd, output_state.buffer, bytes, offset) != 0)
            return -1;
#endif
        output_state.used = 0;
        output_state.last_chunk = id;
#else
        if (swbwa_output_flush() != 0) return -1;
#endif
        if (output_state.debug_enabled) {
            swbwa_output_extent_t *extent = &output_state.extents[output_state.extent_count++];
            extent->id = id;
            extent->start = start;
            extent->end = end;
            extent->reads = reads;
            extent->bytes = bytes;
            extent->offset = bytes || SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
                           ? output_state.last_offset : 0;
        }
#if SWBWA_OUTPUT_RMA_ONLY
        if (swbwa_output_end_chunk() != 0) return -1;
#endif
        return 0;
    }
#else
    (void)sam_lengths;
    for (i = 0; i < reads; ++i) {
        if (seqs[i].sam && swbwa_output_write(seqs[i].sam, strlen(seqs[i].sam)) != 0)
            return -1;
    }
    return swbwa_output_flush();
#endif
}

int swbwa_output_flush(void)
{
    size_t bytes;
    double start;
    int result;

    if (!output_state.opened) {
        errno = EINVAL;
        return -1;
    }
    bytes = output_state.used;
    if (bytes == 0) return 0;
    start = output_debug_now();
    if (output_state.debug_enabled) {
        ++output_state.buffered_flush_calls;
        output_state.buffered_flush_bytes += (uint64_t)bytes;
    }
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    errno = EINVAL;
    result = -1;
#elif SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    result = flush_single_unordered();
#elif SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    output_state.used = 0;
    result = 0;
#else
    result = write_all(output_state.fd, output_state.buffer, bytes);
    if (result == 0) output_state.used = 0;
#endif
    if (output_state.debug_enabled)
        output_state.buffered_flush_seconds += output_debug_now() - start;
    return result;
}

static const char *output_debug_mode_name(void)
{
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    return SWBWA_OUTPUT_RMA_ONLY ? "MPI single_ordered RMA-only (no disk writes; sampled validation)"
                                : "MPI single_ordered";
#endif
#if SWBWA_OUTPUT_RMA_ONLY
    return "MPI single_unordered RMA-only (no disk writes; sampled validation)";
#endif
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    return "MPI single_unordered";
#elif SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SPLIT
    return "MPI split";
#elif SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    return SWBWA_USE_MPI ? "MPI discard" : "non-MPI discard";
#else
    return "non-MPI POSIX";
#endif
}

static void print_output_debug_time(const char *label, double seconds,
                                    uint64_t calls, int indent)
{
    fprintf(stderr, "%*s%-42s %10.6f s", indent, "", label, seconds);
    if (calls > 0)
        fprintf(stderr, "  (%" PRIu64 " calls, %9.3f us/call)",
                calls, seconds * 1.0e6 / (double)calls);
    fputc('\n', stderr);
}

static void print_output_debug_rate(const char *label, uint64_t bytes,
                                    double seconds, int indent)
{
    double mib = (double)bytes / (1024.0 * 1024.0);
    double rate = seconds > 0.0 ? mib / seconds : 0.0;

    fprintf(stderr, "%*s%-42s %10.3f MiB/s\n",
            indent, "", label, rate);
}

static void print_output_debug_report_body(void)
{
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    double reservation_overhead =
        output_state.reservation_total_seconds -
        output_state.fetch_and_op_seconds -
        output_state.win_flush_seconds;

    if (reservation_overhead < 0.0) reservation_overhead = 0.0;
#endif

    fprintf(stderr,
            "\n"
            "====================== SWBWA Output Debug =====================\n"
            "  MPI rank: %06d / %06d\n"
            "  mode:     %s\n"
            "  path:     %s\n"
            "\n"
            "  Output volume\n"
            "    output buffer capacity                       %12zu\n"
            "    swbwa_output_write calls                    %12" PRIu64 "\n"
            "    submitted SAM bytes                         %12" PRIu64 "\n"
            "    buffered flushes                            %12" PRIu64 "\n"
            "    buffered flush bytes                        %12" PRIu64 "\n"
            "    direct oversized writes                     %12" PRIu64 "\n"
            "    direct oversized bytes                      %12" PRIu64 "\n",
            swbwa_mpi_rank(), swbwa_mpi_size(),
            output_debug_mode_name(),
            output_state.name != NULL ? output_state.name : "(none)",
            output_state.capacity,
            output_state.write_calls,
            output_state.submitted_bytes,
            output_state.buffered_flush_calls,
            output_state.buffered_flush_bytes,
            output_state.direct_write_calls,
            output_state.direct_write_bytes);

    fprintf(stderr, "\n  Buffering\n");
    print_output_debug_time(
        "buffered flush total",
        output_state.buffered_flush_seconds,
        output_state.buffered_flush_calls, 4);
    print_output_debug_time("chunk length pass", output_state.chunk_measure_seconds,
                            output_state.chunk_calls, 4);
    print_output_debug_time("chunk packing", output_state.chunk_pack_seconds,
                            output_state.chunk_calls, 4);

#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    fprintf(stderr,
            "\n"
            "  Global offset reservation\n"
            "    target rank                                  %12d\n"
            "    target relationship                          %12s\n"
            "    reservations                                 %12" PRIu64 "\n",
            0, swbwa_mpi_rank() == 0 ? "local" : "remote",
            output_state.reservation_calls);
    print_output_debug_time(
        "reservation total",
        output_state.reservation_total_seconds,
        output_state.reservation_calls, 4);
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    fprintf(stderr, "    ordered prefix polls                         %12" PRIu64 "\n"
                    "    ordered empty chunks                         %12" PRIu64 "\n",
            output_state.ordered_polls, output_state.ordered_empty_chunks);
    print_output_debug_time("ordered prefix wait and publication",
                            output_state.ordered_wait_seconds,
                            output_state.reservation_calls, 6);
#else
    print_output_debug_time(
        "MPI_Fetch_and_op",
        output_state.fetch_and_op_seconds,
        output_state.reservation_calls, 6);
    print_output_debug_time(
        "MPI_Win_flush(target rank 0)",
        output_state.win_flush_seconds,
        output_state.reservation_calls, 6);
    print_output_debug_time(
        "reservation call overhead (derived)",
        reservation_overhead,
        output_state.reservation_calls, 6);
#endif

    fprintf(stderr,
            "\n"
            "  POSIX positioned writes\n"
            "    pwrite system calls                          %12" PRIu64 "\n"
            "    pwrite system-call bytes                     %12" PRIu64 "\n",
            output_state.pwrite_calls,
            output_state.pwrite_bytes);
    print_output_debug_time(
        "POSIX pwrite system calls",
        output_state.pwrite_seconds,
        output_state.pwrite_calls, 4);
    print_output_debug_rate(
        "POSIX pwrite effective bandwidth",
        output_state.pwrite_bytes,
        output_state.pwrite_seconds, 4);
#else
    fprintf(stderr,
            "\n"
            "  POSIX file writes\n"
            "    write system calls                           %12" PRIu64 "\n"
            "    write system-call bytes                      %12" PRIu64 "\n",
            output_state.posix_write_calls,
            output_state.posix_write_bytes);
    print_output_debug_time(
        "POSIX write system calls",
        output_state.posix_write_seconds,
        output_state.posix_write_calls, 4);
    print_output_debug_rate(
        "POSIX write effective bandwidth",
        output_state.posix_write_bytes,
        output_state.posix_write_seconds, 4);
#endif

    fprintf(stderr, "\n  Output close timing\n");
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    print_output_debug_time(
        "MPI_Win_unlock_all", output_state.win_unlock_all_seconds, 1, 4);
    print_output_debug_time(
        "MPI_Win_free", output_state.win_free_seconds, 1, 4);
    print_output_debug_time(
        "POSIX close", output_state.file_close_seconds,
        output_state.owns_fd ? 1 : 0, 4);
#else
    print_output_debug_time(
        "POSIX close", output_state.file_close_seconds,
        output_state.owns_fd ? 1 : 0, 4);
#endif
    fprintf(stderr,
            "================================================================\n"
            "\n");
    for (size_t i = 0; i < output_state.extent_count; ++i) {
        const swbwa_output_extent_t *extent = &output_state.extents[i];
        fprintf(stderr,
                "[SWBWA output extent rank %06d/%06d] chunk=%" PRId64
                " start=%" PRId64 " end=%" PRId64 " reads=%d"
                " offset=%" PRIu64 " bytes=%" PRIu64 "\n",
                swbwa_mpi_rank(), swbwa_mpi_size(), extent->id,
                extent->start, extent->end, extent->reads, extent->offset, extent->bytes);
    }
}

static void output_debug_report(void)
{
    if (!output_state.debug_enabled) return;
    swbwa_mpi_print_rank_ordered(print_output_debug_report_body);
}

int swbwa_output_close(void)
{
    int status = 0;
    int result;
    double start;

    if (!output_state.opened) return 0;
#if SWBWA_OUTPUT_RMA_ONLY
    if (output_state.chunk_active) {
        errno = EINVAL;
        return -1;
    }
#endif
    if (swbwa_output_flush() != 0) status = -1;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    {
        uint64_t local = output_state.reservation_calls + output_state.ordered_empty_chunks;
        uint64_t total = 0;
        if (mpi_check(MPI_Allreduce(&local, &total, 1, MPI_UINT64_T, MPI_SUM,
                                    MPI_COMM_WORLD), "validate ordered chunk count") != 0 ||
            total != (uint64_t)output_state.ordered_count) {
            errno = EINVAL;
            status = -1;
        }
    }
#endif
    start = output_debug_now();
    result = MPI_Win_unlock_all(output_state.offset_window);
    if (output_state.debug_enabled)
        output_state.win_unlock_all_seconds += output_debug_now() - start;
    if (mpi_check(result, "MPI_Win_unlock_all") != 0)
        status = -1;

    start = output_debug_now();
    result = MPI_Win_free(&output_state.offset_window);
    if (output_state.debug_enabled)
        output_state.win_free_seconds += output_debug_now() - start;
    if (mpi_check(result, "MPI_Win_free") != 0)
        status = -1;

    if (output_state.owns_fd) {
        start = output_debug_now();
        result = close(output_state.fd);
        if (output_state.debug_enabled)
            output_state.file_close_seconds += output_debug_now() - start;
        if (result != 0) status = -1;
    }
#else
    if (output_state.owns_fd) {
        start = output_debug_now();
        result = close(output_state.fd);
        if (output_state.debug_enabled)
            output_state.file_close_seconds += output_debug_now() - start;
        if (result != 0) status = -1;
    }
#endif

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    print_discard_hash();
#else
    output_debug_report();
#endif
    if (output_state.stream_md5) {
        if (status == 0) status = swbwa_sam_md5_close(stderr);
        if (status != 0) swbwa_sam_md5_abort();
    }
#if SWBWA_OUTPUT_RMA_ONLY
    swbwa_mpi_print_rank_ordered(print_output_samples);
    free(output_state.samples);
#endif
    free(output_state.name);
    free(output_state.buffer);
    free(output_state.record_lengths);
    free(output_state.extents);
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    return status;
}

const char *swbwa_output_name(void)
{
    return output_state.name;
}
