#include "swbwa_config.h"
#include "swbwa_mpi.h"
#include "swbwa_output.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
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

typedef struct {
    unsigned char *buffer;
    size_t used;
    size_t capacity;
    uint64_t write_calls;
    uint64_t submitted_bytes;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
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
    uint64_t pwrite_calls;
    uint64_t pwrite_bytes;
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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    int discard_hash_enabled;
#endif
    char *name;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    MPI_Win offset_window;
    uint64_t *offset_base;
#endif
} swbwa_output_state_t;

static swbwa_output_state_t output_state;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
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

static uint64_t hash_sam_record(const void *data, size_t length)
{
    static const uint64_t multiplier = UINT64_C(0xc6a4a7935bd1e995);
    static const unsigned int shift = 47;
    const unsigned char *cursor = data;
    size_t remaining = length;
    uint64_t hash = UINT64_C(0x9e3779b97f4a7c15) ^
                    ((uint64_t)length * multiplier);

    while (remaining >= sizeof(uint64_t)) {
        uint64_t word;

        memcpy(&word, cursor, sizeof(word));
        word *= multiplier;
        word ^= word >> shift;
        word *= multiplier;
        hash ^= word;
        hash *= multiplier;
        cursor += sizeof(word);
        remaining -= sizeof(word);
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
}

static void print_discard_hash(void)
{
    fprintf(stderr,
            "[SWBWA output hash rank %06d/%06d] calls=%" PRIu64
            " bytes=%" PRIu64 " sum=0x%016" PRIx64
            " xor=0x%016" PRIx64 " enabled=%d\n",
            swbwa_mpi_rank(), swbwa_mpi_size(),
            output_state.write_calls, output_state.submitted_bytes,
            output_state.hash_sum, output_state.hash_xor,
            output_state.discard_hash_enabled);
    fflush(stderr);
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

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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

    /*
     * RMA gives every rank a disjoint file extent.  POSIX pwrite keeps those
     * writes parallel without routing bulk I/O through Sunway MPI progress.
     */
    return pwrite_all(output_state.fd, data, length, offset);
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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    size_t capacity = 0;
#else
    size_t capacity = (size_t)SWBWA_OUTPUT_BUFFER_BYTES;
#endif

    if (output_state.opened) {
        errno = EALREADY;
        return -1;
    }
#if !(SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)
    if (capacity == 0 || capacity > INT_MAX) {
        errno = EINVAL;
        return -1;
    }
#endif
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    output_state.capacity = capacity;
    output_state.debug_enabled = debug_enabled != 0;
#if !(SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)
    output_state.buffer = malloc(capacity);
    if (output_state.buffer == NULL) return -1;
#endif

#if SWBWA_USE_MPI
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    (void)path;
    output_state.discard_hash_enabled = discard_hash_enabled();
    if (output_state.discard_hash_enabled < 0) goto fail;
    output_state.name = strdup("(discard)");
    if (output_state.name == NULL) goto fail;
#else
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
#elif SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    int local_open_failed;
    int local_open_errno;
    int any_open_failed;
    int open_flags = O_CREAT | O_WRONLY;

    if (sizeof(off_t) < sizeof(int64_t)) {
        errno = EOVERFLOW;
        goto fail;
    }
    output_state.name = strdup(path);
    if (output_state.name == NULL) goto fail;
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
    if (mpi_check(MPI_Win_allocate(swbwa_mpi_rank() == 0 ? sizeof(uint64_t) : 0,
                                   sizeof(uint64_t), MPI_INFO_NULL, MPI_COMM_WORLD,
                                   &output_state.offset_base,
                                   &output_state.offset_window),
                  "MPI_Win_allocate") != 0)
        goto fail;
    if (mpi_check(MPI_Win_lock_all(0, output_state.offset_window),
                  "MPI_Win_lock_all") != 0)
        goto fail_window;
    if (swbwa_mpi_rank() == 0) *output_state.offset_base = 0;
    if (mpi_check(MPI_Win_sync(output_state.offset_window), "MPI_Win_sync") != 0)
        goto fail_locked_window;
    if (mpi_check(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier") != 0)
        goto fail_locked_window;
#endif
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

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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

int swbwa_output_write(const void *data, size_t length)
{
    const unsigned char *source = data;

    if (!output_state.opened || (data == NULL && length != 0)) {
        errno = EINVAL;
        return -1;
    }
    if (length == 0) return 0;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    {
        ++output_state.write_calls;
        output_state.submitted_bytes += (uint64_t)length;
        if (output_state.discard_hash_enabled) {
            uint64_t hash = hash_sam_record(source, length);

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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    result = flush_single_unordered();
#elif SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    return "MPI single_unordered";
#elif SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SPLIT
    return "MPI split";
#elif SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    return "MPI discard";
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

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    print_output_debug_time(
        "MPI_Win_unlock_all", output_state.win_unlock_all_seconds, 1, 4);
    print_output_debug_time(
        "MPI_Win_free", output_state.win_free_seconds, 1, 4);
    print_output_debug_time(
        "POSIX close", output_state.file_close_seconds, 1, 4);
#else
    print_output_debug_time(
        "POSIX close", output_state.file_close_seconds,
        output_state.owns_fd ? 1 : 0, 4);
#endif
    fprintf(stderr,
            "================================================================\n"
            "\n");
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
    if (swbwa_output_flush() != 0) status = -1;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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

    start = output_debug_now();
    result = close(output_state.fd);
    if (output_state.debug_enabled)
        output_state.file_close_seconds += output_debug_now() - start;
    if (result != 0) status = -1;
#else
    if (output_state.owns_fd) {
        start = output_debug_now();
        result = close(output_state.fd);
        if (output_state.debug_enabled)
            output_state.file_close_seconds += output_debug_now() - start;
        if (result != 0) status = -1;
    }
#endif

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD
    print_discard_hash();
#else
    output_debug_report();
#endif
    free(output_state.name);
    free(output_state.buffer);
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    return status;
}

const char *swbwa_output_name(void)
{
    return output_state.name;
}
