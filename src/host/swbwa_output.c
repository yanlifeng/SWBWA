#include "swbwa_config.h"
#include "swbwa_mpi.h"
#include "swbwa_output.h"
#include "swbwa_discard_digest.h"
#include "swbwa_host_workers.h"
#include "swbwa_sam_md5.h"
#if SWBWA_ENABLE_GZIP_OUTPUT
#include "swbwa_gzip_output.h"
#include <sys/uio.h>
#endif

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
#include <athread.h>
#endif
#if SWBWA_OUTPUT_RMA_ONLY
#include <zlib.h>
#endif

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

#if SWBWA_ENABLE_GZIP_OUTPUT
enum { SWBWA_GZIP_IOV_CAPACITY = 1024 };
#endif

#if SWBWA_OUTPUT_RMA_ONLY
enum { OUTPUT_SAMPLE_READS = 100 };
typedef struct {
    int64_t id, start, end;
    int expected_reads, paired;
    uint64_t reads, bytes, sample_reads, sample_bytes, sam_records, sum, xor;
    uint32_t packed_crc32;
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    swbwa_gzip_output *gzip;
    uint64_t gzip_batches;
    uint64_t gzip_file_bytes;
    int gzip_parallel_write;
    int gzip_prepare_output;
    int gzip_prefix_ahead;
    unsigned gzip_iov_limit;
    double gzip_write_wall_seconds;
#endif
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    int discard_hash_enabled;
#endif
    char *name;
#if SWBWA_OUTPUT_RMA_ONLY
    int sample_limit;
    int chunk_active;
    swbwa_output_sample_t chunk;
    swbwa_output_sample_t *samples;
    size_t sample_count, sample_capacity;
#endif
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    MPI_Win offset_window;
    uint64_t *offset_base;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    int ordered_online;
    int64_t ordered_count, last_chunk;
    uint64_t ordered_empty_chunks, ordered_polls;
    unsigned ordered_poll_us;
    double ordered_wait_seconds;
#endif
#endif
} swbwa_output_state_t;

static swbwa_output_state_t output_state;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED && SWBWA_ENABLE_GZIP_OUTPUT
typedef struct prefix_request {
    void *handle;
    int64_t id;
    uint64_t bytes, offset;
    int ready;
    struct prefix_request *next;
} prefix_request_t;

static struct {
    pthread_mutex_t mutex;
    pthread_cond_t changed;
    pthread_t thread;
    prefix_request_t requests[SWBWA_PIPELINE_BUFFER_COUNT], *first, *last;
    int started, initialized, stopping, error;
    int64_t last_id;
    uint64_t submitted, idle_probes;
    double writer_wait_s, idle_probe_s;
} prefix_ahead = {.mutex = PTHREAD_MUTEX_INITIALIZER, .changed = PTHREAD_COND_INITIALIZER};

static void prefix_check(int error)
{
    if (!error) return;
    fprintf(stderr, "[E::output prefix] synchronization failed: %s\n", strerror(error));
    abort();
}

static prefix_request_t *prefix_find(void *handle)
{
    unsigned i;
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        if (prefix_ahead.requests[i].handle == handle) return &prefix_ahead.requests[i];
    return NULL;
}
#endif

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
                " sample_reads=%" PRIu64 " sample_bytes=%" PRIu64 " sam_records=%" PRIu64
                " sum=0x%016" PRIx64 " xor=0x%016" PRIx64
                " packed_crc32=0x%08" PRIx32 "\n",
                swbwa_mpi_rank(), swbwa_mpi_size(), s->id, s->start, s->end,
                s->paired, s->reads, s->bytes, s->sample_reads, s->sample_bytes, s->sam_records,
                s->sum, s->xor, s->packed_crc32);
    }
    fprintf(stderr,
            "[SWBWA output RMA-only rank %06d/%06d] chunks=%zu"
            " reads=%" PRIu64 " sam_bytes=%" PRIu64
            " reserved_bytes=%" PRIu64 " reservations=%" PRIu64
            " disk_bytes=%" PRIu64 " sample_limit=%d\n",
            swbwa_mpi_rank(), swbwa_mpi_size(), output_state.sample_count,
            output_state.write_calls, output_state.submitted_bytes,
            output_state.reserved_bytes, output_state.reservation_calls,
            output_state.pwrite_bytes, output_state.sample_limit);
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
    int gzip = length >= 3 && strcmp(path + length - 3, ".gz") == 0;

    if (gzip) stem_length -= 3;

    if (stem_length >= sizeof(suffix) - 1 &&
        memcmp(path + stem_length - (sizeof(suffix) - 1), suffix, sizeof(suffix) - 1) == 0)
        stem_length -= sizeof(suffix) - 1;
    capacity = stem_length + sizeof(".rank000000.sam.gz") + 16;
    name = malloc(capacity);
    if (name == NULL) return NULL;
    snprintf(name, capacity, "%.*s.rank%06d.sam%s", (int)stem_length, path, rank, gzip ? ".gz" : "");
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
static void ordered_poll_pause(uint64_t polls)
{
    struct timespec deadline;
    if (polls <= 8 || output_state.ordered_poll_us == 0 ||
        clock_gettime(CLOCK_MONOTONIC, &deadline) != 0) {
        sched_yield();
        return;
    }
    deadline.tv_nsec += (long)output_state.ordered_poll_us * 1000L;
    if (deadline.tv_nsec >= 1000000000L) {
        ++deadline.tv_sec;
        deadline.tv_nsec -= 1000000000L;
    }
    /* Absolute deadlines survive the frequent signals on Sunway. No local
     * dispatch/runtime mutex is held while this writer yields its MPE. */
    while (clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &deadline, NULL) == EINTR) {}
}

/* Global input tickets make each rank's FIFO monotonic. Therefore the first
 * unpublished prefix cannot be blocked behind a later chunk on that rank. */
static int reserve_ordered(int64_t id, uint64_t bytes, uint64_t *offset,
                           int record_stats)
{
    uint64_t next, polls = 0;
    double start = output_debug_now();
    int result;

    if (id < 0 || id == INT64_MAX ||
        (!output_state.ordered_online && id >= output_state.ordered_count) || bytes > INT64_MAX) {
        errno = EINVAL;
        return -1;
    }
    if (output_state.ordered_online) {
        uint64_t ready;
        /* Exactly one owner publishes each ID. Flush the new byte offset
         * before its ID, so the next owner cannot observe a stale prefix. */
        do {
            result = MPI_Fetch_and_op(NULL, &ready, MPI_UINT64_T, 0, 0,
                                      MPI_NO_OP, output_state.offset_window);
            if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
            if (mpi_check(result, "read online ordered ID") != 0) return -1;
            ++polls;
            if (ready > (uint64_t)id) { errno = EINVAL; return -1; }
            if (ready != (uint64_t)id) ordered_poll_pause(polls);
        } while (ready != (uint64_t)id);
        result = MPI_Fetch_and_op(NULL, offset, MPI_UINT64_T, 0, 1,
                                  MPI_NO_OP, output_state.offset_window);
        if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
        if (mpi_check(result, "read online ordered offset") != 0) return -1;
    } else do {
        result = MPI_Fetch_and_op(NULL, offset, MPI_UINT64_T, 0, (MPI_Aint)id,
                                  MPI_NO_OP, output_state.offset_window);
        if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
        if (mpi_check(result, "read ordered prefix") != 0) return -1;
        ++polls;
        if (*offset == UINT64_MAX) ordered_poll_pause(polls);
    } while (*offset == UINT64_MAX);
    if (*offset > (uint64_t)INT64_MAX - bytes) {
        errno = EOVERFLOW;
        return -1;
    }
    next = *offset + bytes;
    /* Accumulate/NO_OP provides atomic access to each prefix slot; ordinary
     * concurrent Put/Get would not provide this atomicity contract. */
    result = MPI_Accumulate(&next, 1, MPI_UINT64_T, 0,
                            output_state.ordered_online ? 1 : (MPI_Aint)(id + 1),
                            1, MPI_UINT64_T, MPI_REPLACE, output_state.offset_window);
    if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
    if (mpi_check(result, "publish ordered prefix") != 0) return -1;
    if (output_state.ordered_online) {
        uint64_t ready = (uint64_t)id + 1;
        result = MPI_Accumulate(&ready, 1, MPI_UINT64_T, 0, 0, 1, MPI_UINT64_T,
                                MPI_REPLACE, output_state.offset_window);
        if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
        if (mpi_check(result, "publish online ordered ID") != 0) return -1;
    }
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

#if SWBWA_ENABLE_GZIP_OUTPUT
static void *prefix_worker(void *unused)
{
    int binding_error = 0;
    (void)unused;
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
    unsigned long actual;
    SET_TO_SWCG(SWBWA_HOST_MPE_THREADS - 1);
    __asm__ volatile("rcid %0" : "=r"(actual));
    if ((actual & 7) != SWBWA_HOST_MPE_THREADS - 1) binding_error = EINVAL;
#endif
    prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
    prefix_ahead.error = binding_error;
    prefix_ahead.initialized = 1;
    prefix_check(pthread_cond_broadcast(&prefix_ahead.changed));
    while (!prefix_ahead.error) {
        prefix_request_t *request;
        int result, saved;
        while (!prefix_ahead.first && !prefix_ahead.stopping) {
            struct timespec deadline;
            int wait_result, available;
            double begin;
            if (swbwa_mpi_progress_thread_active()) {
                prefix_check(pthread_cond_wait(&prefix_ahead.changed, &prefix_ahead.mutex));
                continue;
            }
            if (clock_gettime(CLOCK_REALTIME, &deadline) != 0) {
                prefix_ahead.error = errno;
                break;
            }
            deadline.tv_nsec += 100000L;
            if (deadline.tv_nsec >= 1000000000L) {
                ++deadline.tv_sec;
                deadline.tv_nsec -= 1000000000L;
            }
            wait_result = pthread_cond_timedwait(&prefix_ahead.changed,
                                                  &prefix_ahead.mutex, &deadline);
            if (wait_result != ETIMEDOUT) { prefix_check(wait_result); continue; }
            /* A software-progress MPI must also serve peers while this rank
             * has no new batch, including when its disk writer is blocked. */
            prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
            begin = output_debug_now();
            result = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                                 &available, MPI_STATUS_IGNORE);
            prefix_ahead.idle_probe_s += output_debug_now() - begin;
            ++prefix_ahead.idle_probes;
            prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
            if (result != MPI_SUCCESS) { prefix_ahead.error = EIO; break; }
        }
        if (prefix_ahead.error) {
            prefix_check(pthread_cond_broadcast(&prefix_ahead.changed));
            break;
        }
        if (!prefix_ahead.first) break;
        request = prefix_ahead.first;
        prefix_ahead.first = request->next;
        if (!prefix_ahead.first) prefix_ahead.last = NULL;
        prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
        result = reserve_ordered(request->id, request->bytes, &request->offset, 1);
        saved = errno;
        prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
        request->ready = 1;
        if (result) prefix_ahead.error = saved ? saved : EIO;
        prefix_check(pthread_cond_broadcast(&prefix_ahead.changed));
    }
    prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
    return NULL;
}

static int prefix_start(void)
{
    int error, level;
    if (!output_state.gzip_prefix_ahead) return 0;
    if (MPI_Query_thread(&level) != MPI_SUCCESS || level < MPI_THREAD_MULTIPLE) {
        errno = ENOTSUP;
        return -1;
    }
    memset(prefix_ahead.requests, 0, sizeof(prefix_ahead.requests));
    prefix_ahead.first = prefix_ahead.last = NULL;
    prefix_ahead.initialized = prefix_ahead.stopping = prefix_ahead.error = 0;
    prefix_ahead.last_id = -1;
    prefix_ahead.submitted = 0;
    prefix_ahead.idle_probes = 0;
    prefix_ahead.writer_wait_s = 0;
    prefix_ahead.idle_probe_s = 0;
    error = pthread_create(&prefix_ahead.thread, NULL, prefix_worker, NULL);
    if (error) { errno = error; return -1; }
    prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
    while (!prefix_ahead.initialized)
        prefix_check(pthread_cond_wait(&prefix_ahead.changed, &prefix_ahead.mutex));
    error = prefix_ahead.error;
    prefix_ahead.started = !error;
    prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
    if (error) { prefix_check(pthread_join(prefix_ahead.thread, NULL)); errno = error; return -1; }
    return 0;
}

static int prefix_wait(void *handle, int64_t id, uint64_t bytes, uint64_t *offset)
{
    prefix_request_t *request;
    int error;
    double begin = output_debug_now();
    prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
    request = prefix_find(handle);
    if (!request || request->id != id || request->bytes != bytes) {
        prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
        errno = EINVAL;
        return -1;
    }
    while (!request->ready && !prefix_ahead.error)
        prefix_check(pthread_cond_wait(&prefix_ahead.changed, &prefix_ahead.mutex));
    error = prefix_ahead.error;
    *offset = request->offset;
    memset(request, 0, sizeof(*request));
    prefix_ahead.writer_wait_s += output_debug_now() - begin;
    prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
    if (error) { errno = error; return -1; }
    return 0;
}

static int prefix_stop(void)
{
    int error;
    if (!prefix_ahead.started) return 0;
    prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
    prefix_ahead.stopping = 1;
    prefix_check(pthread_cond_broadcast(&prefix_ahead.changed));
    prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
    prefix_check(pthread_join(prefix_ahead.thread, NULL));
    prefix_ahead.started = 0;
    error = prefix_ahead.error;
    fprintf(stderr, "[SWBWA gzip prefix rank %06d/%06d] ahead=1 submitted=%" PRIu64
            " writer_wait_s=%.6f error=%d idle_probes=%" PRIu64 " idle_probe_s=%.6f\n",
            swbwa_mpi_rank(), swbwa_mpi_size(), prefix_ahead.submitted,
            prefix_ahead.writer_wait_s, error, prefix_ahead.idle_probes,
            prefix_ahead.idle_probe_s);
    if (error) { errno = error; return -1; }
    return 0;
}
#endif

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

static int reserve_unordered(uint64_t length, uint64_t *reserved_offset)
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
    *reserved_offset = offset;
    return 0;
}

static int write_single_unordered(const unsigned char *data, size_t length)
{
    uint64_t offset;
    if (length == 0) return 0;
    if (reserve_unordered(length, &offset) != 0) return -1;

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

static int output_open(const char *path, int debug_enabled, int online)
{
    const char *md5_option = getenv("SWBWA_OUTPUT_MD5");
    int stream_md5 = md5_option != NULL && strcmp(md5_option, "1") == 0;
    size_t path_length = path ? strlen(path) : 0;
    int gzip = path_length >= 3 && strcmp(path + path_length - 3, ".gz") == 0;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    size_t capacity = 0;
#else
    size_t capacity = (size_t)SWBWA_OUTPUT_BUFFER_BYTES;
#endif

    if (output_state.opened) {
        errno = EALREADY;
        return -1;
    }
    if (gzip && (stream_md5 || SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)) {
        fprintf(stderr, "[E::output] gzip output cannot be combined with discard or streaming MD5\n");
        errno = EINVAL;
        return -1;
    }
#if !SWBWA_ENABLE_GZIP_OUTPUT
    if (gzip) {
        fprintf(stderr, "[E::output] .gz output requires GZIP_OUTPUT=1\n");
        errno = ENOTSUP;
        return -1;
    }
#endif
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (gzip) {
        long iov_limit = sysconf(_SC_IOV_MAX);
        const char *parallel = getenv("SWBWA_GZIP_PARALLEL_WRITE");
        const char *prepare = getenv("SWBWA_GZIP_PREPARE_OUTPUT");
        output_state.gzip_iov_limit = iov_limit > 0
            ? (unsigned)(iov_limit < SWBWA_GZIP_IOV_CAPACITY ? iov_limit : SWBWA_GZIP_IOV_CAPACITY)
            : 16; /* POSIX minimum if the runtime limit is unavailable. */
        if (parallel && strcmp(parallel, "0") && strcmp(parallel, "1")) {
            errno = EINVAL;
            goto fail;
        }
        output_state.gzip_parallel_write = parallel && strcmp(parallel, "0");
        if (prepare && strcmp(prepare, "0") && strcmp(prepare, "1")) {
            errno = EINVAL;
            goto fail;
        }
        output_state.gzip_prepare_output = !prepare || !strcmp(prepare, "1");
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        {
            const char *ahead = getenv("SWBWA_GZIP_PREFIX_AHEAD");
            if (ahead && strcmp(ahead, "0") && strcmp(ahead, "1")) { errno = EINVAL; goto fail; }
            output_state.gzip_prefix_ahead = ahead && !strcmp(ahead, "1") &&
                output_state.gzip_prepare_output;
        }
#endif
        output_state.gzip = swbwa_gzip_output_create(SWBWA_CPE_COUNT);
        if (!output_state.gzip) goto fail;
        capacity = 0;
        output_state.capacity = 0;
    }
#endif
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
    {
        const char *value = getenv("SWBWA_ORDERED_POLL_US");
        char *end;
        unsigned long delay = 0;
        if (value) {
            errno = 0;
            delay = strtoul(value, &end, 10);
            if (errno || !*value || *end || delay > 1000) {
                errno = EINVAL;
                goto fail;
            }
        }
        output_state.ordered_poll_us = (unsigned)delay;
    }
    output_state.ordered_online = online;
#else
    (void)online;
#endif
#if SWBWA_OUTPUT_RMA_ONLY
    {
        const char *value = getenv("SWBWA_OUTPUT_SAMPLE_READS");
        output_state.sample_limit = OUTPUT_SAMPLE_READS;
        if (value && *value) {
            char *end;
            long limit;
            errno = 0;
            limit = strtol(value, &end, 10);
            if (errno || *end || limit < 0 || limit > INT_MAX) {
                errno = EINVAL;
                goto fail;
            }
            output_state.sample_limit = (int)limit;
        }
    }
#endif
#if !(SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)
    if (capacity) {
        output_state.buffer = malloc(capacity);
        if (output_state.buffer == NULL) goto fail;
    }
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
    output_state.ordered_count = -1;
    output_state.last_chunk = -1;
    if (!online) {
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
        output_state.ordered_count = swbwa_mpi_fastq_scheduler_chunk_count();
        if (strcmp(swbwa_mpi_fastq_scheduler_ticket_mode(), "global") != 0) {
            errno = EINVAL;
            goto fail;
        }
#endif
        if (output_state.ordered_count < 0 ||
            (uint64_t)output_state.ordered_count >= (uint64_t)PTRDIFF_MAX / sizeof(uint64_t)) {
            errno = EINVAL;
            goto fail;
        }
    }
    window_bytes = online ? 2 * sizeof(uint64_t)
                          : ((size_t)output_state.ordered_count + 1) * sizeof(uint64_t);
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
        memset(output_state.offset_base, online ? 0 : 0xff, window_bytes);
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

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED && SWBWA_ENABLE_GZIP_OUTPUT
    if (prefix_start() != 0) goto fail_locked_window;
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    swbwa_gzip_output_destroy(output_state.gzip);
#endif
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    return -1;
}

int swbwa_output_open(const char *path, int debug_enabled)
{
    return output_open(path, debug_enabled, 0);
}

#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
int swbwa_output_open_online(const char *path, int debug_enabled)
{
    return output_open(path, debug_enabled, 1);
}

int swbwa_output_ordered_finish(int64_t chunks)
{
    if (!output_state.opened || !output_state.ordered_online || chunks < 0 ||
        output_state.ordered_count != -1) {
        errno = EINVAL;
        return -1;
    }
    output_state.ordered_count = chunks;
    return 0;
}
#endif

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

#if SWBWA_ENABLE_GZIP_OUTPUT
typedef struct {
    const struct iovec *vectors;
    unsigned count, workers;
    unsigned first[SWBWA_HOST_MPE_THREADS], end[SWBWA_HOST_MPE_THREADS];
    uint64_t offsets[SWBWA_HOST_MPE_THREADS];
    uint64_t calls[SWBWA_HOST_MPE_THREADS], bytes[SWBWA_HOST_MPE_THREADS];
    double seconds[SWBWA_HOST_MPE_THREADS];
    int errors[SWBWA_HOST_MPE_THREADS];
} gzip_write_task;

/* Scatter/gather writes preserve member order without packing compressed data.
 * Short writes and EINTR may stop inside a member and are resumed exactly. */
static void write_gzip_slice(void *opaque, int worker, int workers)
{
    gzip_write_task *task = opaque;
    const struct iovec *vectors = task->vectors;
    unsigned first = task->first[worker], end = task->end[worker];
    uint64_t offset = task->offsets[worker];
    size_t skip = 0;
    if (workers == 1 && task->workers > 1) {
        first = 0;
        end = task->count;
        offset = task->offsets[0];
    }
    while (first < end) {
        struct iovec batch[SWBWA_GZIP_IOV_CAPACITY];
        unsigned n = end - first > output_state.gzip_iov_limit
                   ? output_state.gzip_iov_limit : end - first, i;
        ssize_t written;
        size_t remaining;
        double start;
        for (i = 0; i < n; ++i) batch[i] = vectors[first + i];
        batch[0].iov_base = (unsigned char *)batch[0].iov_base + skip;
        batch[0].iov_len -= skip;
        start = output_debug_now();
        written = pwritev(output_state.fd, batch, (int)n, (off_t)offset);
        if (output_state.debug_enabled) {
            ++task->calls[worker];
            task->seconds[worker] += output_debug_now() - start;
            if (written > 0) task->bytes[worker] += (uint64_t)written;
        }
        if (written < 0) { if (errno == EINTR) continue; task->errors[worker] = errno; return; }
        if (!written) { task->errors[worker] = EIO; return; }
        offset += (uint64_t)written;
        remaining = (size_t)written;
        while (first < end && remaining >= vectors[first].iov_len - skip) {
            remaining -= vectors[first].iov_len - skip;
            ++first;
            skip = 0;
        }
        skip += remaining;
    }
}

static int write_gzip_vectors(const struct iovec *vectors, unsigned count, uint64_t offset)
{
#if SWBWA_OUTPUT_RMA_ONLY
    (void)vectors; (void)count; (void)offset;
    return 0;
#else
    gzip_write_task task;
    unsigned worker, index = 0;
    unsigned workers = output_state.gzip_parallel_write && count >= SWBWA_HOST_MPE_THREADS
                     ? SWBWA_HOST_MPE_THREADS : 1;
    double start = output_debug_now();
    memset(&task, 0, sizeof(task));
    task.vectors = vectors;
    task.count = count;
    task.workers = workers;
    for (worker = 0; worker < workers; ++worker) {
        task.first[worker] = index;
        task.end[worker] = (unsigned)((uint64_t)count * (worker + 1) / workers);
        task.offsets[worker] = offset;
        for (; index < task.end[worker]; ++index) offset += vectors[index].iov_len;
    }
    if (workers > 1) {
        if (swbwa_host_io_workers_run(write_gzip_slice, &task) != 0) return -1;
    } else write_gzip_slice(&task, 0, 1);
    output_state.gzip_write_wall_seconds += output_debug_now() - start;
    for (worker = 0; worker < workers; ++worker) {
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
        output_state.pwrite_calls += task.calls[worker];
        output_state.pwrite_bytes += task.bytes[worker];
        output_state.pwrite_seconds += task.seconds[worker];
#else
        output_state.posix_write_calls += task.calls[worker];
        output_state.posix_write_bytes += task.bytes[worker];
        output_state.posix_write_seconds += task.seconds[worker];
#endif
        if (task.errors[worker]) { errno = task.errors[worker]; return -1; }
    }
    return 0;
#endif
}

static int gzip_stream_write(const void *data, size_t length)
{
    const struct iovec *vectors;
    unsigned count;
    uint64_t bytes, offset = output_state.gzip_file_bytes;
    if (swbwa_gzip_output_bytes(output_state.gzip, data, length, &vectors, &count, &bytes) != 0)
        return -1;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    if (reserve_unordered(bytes, &offset) != 0) return -1;
#endif
    if (bytes > INT64_MAX || offset > (uint64_t)INT64_MAX - bytes) {
        errno = EOVERFLOW;
        return -1;
    }
    if (write_gzip_vectors(vectors, count, offset) != 0) return -1;
    ++output_state.gzip_batches;
    output_state.gzip_file_bytes += bytes;
    return 0;
}

static int gzip_record_extent(int64_t id, int64_t start, int64_t end,
                              int reads, uint64_t offset, uint64_t bytes)
{
    swbwa_output_extent_t *extent;
    if (!output_state.debug_enabled) return 0;
    if (output_state.extent_count == output_state.extent_capacity) {
        size_t capacity = output_state.extent_capacity ? output_state.extent_capacity * 2 : 16;
        swbwa_output_extent_t *extents;
        if (capacity < output_state.extent_capacity || capacity > SIZE_MAX / sizeof(*extents)) {
            errno = EOVERFLOW;
            return -1;
        }
        extents = realloc(output_state.extents, capacity * sizeof(*extents));
        if (!extents) return -1;
        output_state.extents = extents;
        output_state.extent_capacity = capacity;
    }
    extent = &output_state.extents[output_state.extent_count++];
    extent->id = id;
    extent->start = start;
    extent->end = end;
    extent->reads = reads;
    extent->offset = offset;
    extent->bytes = bytes;
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
        if (!output_state.sample_limit || s->reads < (uint64_t)output_state.sample_limit) {
            uint64_t hash = hash_sam_record(data, length);
            const unsigned char *cursor = data, *limit = cursor + length;
            while (cursor < limit) {
                const unsigned char *lf = memchr(cursor, '\n', (size_t)(limit - cursor));
                if (!lf) break;
                ++s->sam_records;
                cursor = lf + 1;
            }
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip) return gzip_stream_write(data, length);
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

int swbwa_output_prepare_chunk(const bseq1_t *seqs, const int *sam_lengths,
                               int reads, void **handle)
{
    if (!handle) { errno = EINVAL; return -1; }
    *handle = NULL;
    if (!output_state.opened || reads < 0 || (reads && !seqs)) {
        errno = EINVAL;
        return -1;
    }
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip && output_state.gzip_prepare_output) {
        swbwa_gzip_output_batch *batch;
        if (swbwa_gzip_output_prepare(output_state.gzip, seqs, sam_lengths, reads, &batch))
            return -1;
        *handle = batch;
    }
#else
    (void)sam_lengths;
#endif
    return 0;
}

int swbwa_output_release_prepared_chunk(void *handle)
{
    if (!handle) return 0;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED && SWBWA_ENABLE_GZIP_OUTPUT
    if (prefix_ahead.started) {
        int pending;
        prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
        pending = prefix_find(handle) != NULL;
        prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
        if (pending) { errno = EBUSY; return -1; }
    }
#endif
#if SWBWA_ENABLE_GZIP_OUTPUT
    return swbwa_gzip_output_release(output_state.gzip, handle);
#else
    errno = EINVAL;
    return -1;
#endif
}

int swbwa_output_publish_prepared_chunk(int64_t id, int reads, void *handle)
{
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED && SWBWA_ENABLE_GZIP_OUTPUT
    if (prefix_ahead.started && handle) {
        prefix_request_t *request;
        const struct iovec *vectors;
        unsigned count;
        uint64_t raw_bytes, bytes;
        int error;
        if (id < 0 || swbwa_gzip_output_batch_data(output_state.gzip, handle, reads,
                &vectors, &count, &raw_bytes, &bytes)) { errno = EINVAL; return -1; }
        prefix_check(pthread_mutex_lock(&prefix_ahead.mutex));
        error = prefix_ahead.error;
        request = prefix_find(NULL);
        if (!error && (id <= prefix_ahead.last_id || prefix_find(handle))) error = EINVAL;
        if (!error && !request) error = EAGAIN;
        if (!error) {
            *request = (prefix_request_t){.handle=handle, .id=id, .bytes=bytes};
            if (prefix_ahead.last) prefix_ahead.last->next = request;
            else prefix_ahead.first = request;
            prefix_ahead.last = request;
            prefix_ahead.last_id = id;
            ++prefix_ahead.submitted;
            prefix_check(pthread_cond_broadcast(&prefix_ahead.changed));
        }
        prefix_check(pthread_mutex_unlock(&prefix_ahead.mutex));
        if (error) { errno = error; return -1; }
    }
#else
    (void)id; (void)reads; (void)handle;
#endif
    return 0;
}

static int output_write_chunk(int64_t id, int64_t start, int64_t end,
                               const bseq1_t *seqs, const int *sam_lengths,
                               int reads, int paired, void *prepared)
{
    int i;
    if (!output_state.opened || id < 0 || start < 0 || end < start ||
        reads < 0 || (reads && seqs == NULL && (!prepared || SWBWA_OUTPUT_RMA_ONLY)) ||
        (paired != 0 && paired != 1) || (paired && (reads & 1))) {
        errno = EINVAL;
        return -1;
    }
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip) {
        const struct iovec *vectors;
        unsigned count;
        uint64_t raw_bytes, gzip_bytes, offset = output_state.gzip_file_bytes;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (id <= output_state.last_chunk) { errno = EINVAL; return -1; }
#endif
        if (prepared) {
            if (swbwa_gzip_output_batch_data(output_state.gzip, prepared, reads,
                                             &vectors, &count, &raw_bytes, &gzip_bytes))
                return -1;
        } else if (swbwa_gzip_output_encode(output_state.gzip, seqs, sam_lengths, reads,
                                             &vectors, &count, &raw_bytes, &gzip_bytes))
            return -1;
#if SWBWA_OUTPUT_RMA_ONLY
        if (swbwa_output_begin_chunk(id, start, end, reads, paired) != 0) return -1;
        for (i = 0; i < reads; ++i) {
            size_t length = sam_lengths ? (size_t)sam_lengths[i] : strlen(seqs[i].sam);
            if (!output_state.sample_limit || i < output_state.sample_limit) {
                const unsigned char *cursor = (const unsigned char *)seqs[i].sam;
                size_t remaining = length;
                uint64_t hash = hash_sam_record(cursor, length);
                while (remaining) {
                    const unsigned char *lf = memchr(cursor, '\n', remaining);
                    if (!lf) break;
                    ++output_state.chunk.sam_records;
                    remaining -= (size_t)(lf - cursor) + 1;
                    cursor = lf + 1;
                }
                ++output_state.chunk.sample_reads;
                output_state.chunk.sample_bytes += length;
                output_state.chunk.sum += hash;
                output_state.chunk.xor ^= hash;
            }
            if (!output_state.sample_limit) {
                size_t pos = 0;
                while (pos < length) {
                    uInt part = length - pos > UINT_MAX ? UINT_MAX : (uInt)(length - pos);
                    output_state.chunk.packed_crc32 = (uint32_t)crc32(output_state.chunk.packed_crc32,
                                            (const Bytef *)seqs[i].sam + pos, part);
                    pos += part;
                }
            }
        }
        output_state.chunk.reads = (uint64_t)reads;
        output_state.chunk.bytes = raw_bytes;
#endif
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (prefix_ahead.started) {
            if (!prepared) { errno = EINVAL; return -1; }
            if (prefix_wait(prepared, id, gzip_bytes, &offset) != 0) return -1;
        } else if (reserve_ordered(id, gzip_bytes, &offset, 1) != 0) return -1;
        output_state.last_chunk = id;
#else
        if (reserve_unordered(gzip_bytes, &offset) != 0) return -1;
#endif
#endif
        if (gzip_bytes > INT64_MAX || offset > (uint64_t)INT64_MAX - gzip_bytes) {
            errno = EOVERFLOW;
            return -1;
        }
        if (write_gzip_vectors(vectors, count, offset) != 0) return -1;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
        if (gzip_record_extent(id, start, end, reads, offset, gzip_bytes) != 0) return -1;
#endif
        ++output_state.gzip_batches;
        output_state.gzip_file_bytes += gzip_bytes;
        ++output_state.chunk_calls;
        output_state.write_calls += (uint64_t)reads;
        output_state.submitted_bytes += raw_bytes;
#if SWBWA_OUTPUT_RMA_ONLY
        if (swbwa_output_end_chunk() != 0) return -1;
#endif
        return 0;
    }
#endif
    if (prepared) { errno = EINVAL; return -1; }
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
    {
        output_pack_t task;
        size_t bytes = 0;
        double clock_start;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (id <= output_state.last_chunk) {
            errno = EINVAL;
            return -1;
        }
#endif

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
        size_t packed_pos = 0;
        for (i = 0; i < reads && (!output_state.sample_limit || i < output_state.sample_limit); ++i) {
            const char *cursor = (const char *)task.packed + packed_pos;
            const char *limit = cursor + task.lengths[i];
            uint64_t hash = hash_sam_record(cursor, task.lengths[i]);
            while (cursor < limit) {
                const char *lf = memchr(cursor, '\n', (size_t)(limit - cursor));
                if (!lf) break;
                ++output_state.chunk.sam_records;
                cursor = lf + 1;
            }
            ++output_state.chunk.sample_reads;
            output_state.chunk.sample_bytes += task.lengths[i];
            output_state.chunk.sum += hash;
            output_state.chunk.xor ^= hash;
            packed_pos += task.lengths[i];
        }
        if (!output_state.sample_limit) {
            size_t pos = 0;
            uLong crc = crc32(0, Z_NULL, 0);
            while (pos < bytes) {
                uInt part = bytes - pos > UINT_MAX ? UINT_MAX : (uInt)(bytes - pos);
                crc = crc32(crc, task.packed + pos, part);
                pos += part;
            }
            output_state.chunk.packed_crc32 = (uint32_t)crc;
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

int swbwa_output_write_prepared_chunk(int64_t id, int64_t start, int64_t end,
                                      const bseq1_t *seqs, const int *sam_lengths,
                                      int reads, int paired, void *handle)
{
    int result = output_write_chunk(id, start, end, seqs, sam_lengths, reads, paired, handle);
    int saved = errno;
    int released = swbwa_output_release_prepared_chunk(handle);
    if (result) errno = saved;
    return result ? result : released;
}

int swbwa_output_write_chunk(int64_t id, int64_t start, int64_t end,
                             const bseq1_t *seqs, const int *sam_lengths,
                             int reads, int paired)
{
    return swbwa_output_write_prepared_chunk(id, start, end, seqs, sam_lengths, reads, paired, NULL);
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
    return SWBWA_OUTPUT_RMA_ONLY ? "MPI single_ordered RMA-only (no disk writes; configurable validation)"
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
                    "    ordered empty chunks                         %12" PRIu64 "\n"
                    "    ordered poll pause after 8 misses (us)       %12u\n",
            output_state.ordered_polls, output_state.ordered_empty_chunks,
            output_state.ordered_poll_us);
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
        "pwrite bytes / summed syscall time",
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
        "write bytes / summed syscall time",
        output_state.posix_write_bytes,
        output_state.posix_write_seconds, 4);
#endif
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip) {
        fprintf(stderr, "    Parallel syscall times are worker sums, not wall time.\n");
        print_output_debug_rate("gzip write wall bandwidth",
            output_state.pwrite_bytes + output_state.posix_write_bytes,
            output_state.gzip_write_wall_seconds, 4);
    }
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip && swbwa_gzip_output_pending(output_state.gzip)) {
        errno = EBUSY;
        return -1;
    }
#endif
#if SWBWA_OUTPUT_RMA_ONLY
    if (output_state.chunk_active) {
        errno = EINVAL;
        return -1;
    }
#endif
    if (swbwa_output_flush() != 0) status = -1;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED && SWBWA_ENABLE_GZIP_OUTPUT
    if (prefix_stop() != 0) status = -1;
#endif
#if SWBWA_ENABLE_GZIP_OUTPUT
    if (output_state.gzip) {
        uint64_t batches = output_state.gzip_batches;
        static const unsigned char empty_gzip[20] = {
            31, 139, 8, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
#if SWBWA_USE_MPI && SWBWA_OUTPUT_SINGLE_FILE
        uint64_t total = 0;
        if (mpi_check(MPI_Allreduce(&batches, &total, 1, MPI_UINT64_T, MPI_SUM,
                                    MPI_COMM_WORLD), "count gzip output batches") != 0) status = -1;
        batches = swbwa_mpi_rank() == 0 ? total : 1;
#endif
#if !SWBWA_OUTPUT_RMA_ONLY
        if (!batches && status == 0 && write_all(output_state.fd, empty_gzip, sizeof(empty_gzip)) != 0)
            status = -1;
#else
        (void)batches;
        (void)empty_gzip;
#endif
        swbwa_gzip_output_report(output_state.gzip, stderr, swbwa_mpi_rank(), swbwa_mpi_size());
        fprintf(stderr, "[SWBWA gzip sink rank %06d/%06d] parallel_write=%d iov_limit=%u write_wall_s=%.6f"
                " syscall_worker_s=%.6f\n", swbwa_mpi_rank(), swbwa_mpi_size(),
                output_state.gzip_parallel_write, output_state.gzip_iov_limit, output_state.gzip_write_wall_seconds,
                output_state.posix_write_seconds + output_state.pwrite_seconds);
    }
#endif

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
        if (output_state.ordered_online) {
            uint64_t local_bytes = output_state.reserved_bytes, total_bytes = 0;
            uint64_t ready = 0, bytes = 0;
            int invalid, any_invalid;
            result = MPI_Allreduce(&local_bytes, &total_bytes, 1, MPI_UINT64_T,
                                    MPI_SUM, MPI_COMM_WORLD);
            if (mpi_check(result, "validate online ordered bytes") != 0) status = -1;
            result = MPI_Fetch_and_op(NULL, &ready, MPI_UINT64_T, 0, 0,
                                      MPI_NO_OP, output_state.offset_window);
            if (result == MPI_SUCCESS)
                result = MPI_Fetch_and_op(NULL, &bytes, MPI_UINT64_T, 0, 1,
                                          MPI_NO_OP, output_state.offset_window);
            if (result == MPI_SUCCESS) result = MPI_Win_flush(0, output_state.offset_window);
            if (mpi_check(result, "validate online ordered final prefix") != 0 ||
                ready != total || bytes != total_bytes) {
                errno = EINVAL;
                status = -1;
            }
            invalid = status != 0;
            result = MPI_Allreduce(&invalid, &any_invalid, 1, MPI_INT, MPI_MAX,
                                    MPI_COMM_WORLD);
            if (mpi_check(result, "validate online ordered status") != 0 || any_invalid) {
                errno = EINVAL;
                status = -1;
            }
            fprintf(stderr, "[SWBWA ordered stream rank %06d/%06d] chunks=%" PRIu64
                    " expected_chunks=%" PRId64 " file_bytes=%" PRIu64 " validated=%d\n",
                    swbwa_mpi_rank(), swbwa_mpi_size(), total,
                    output_state.ordered_count, bytes, status == 0);
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
#if SWBWA_ENABLE_GZIP_OUTPUT
    swbwa_gzip_output_destroy(output_state.gzip);
#endif
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    return status;
}

const char *swbwa_output_name(void)
{
    return output_state.name;
}
