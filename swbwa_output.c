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
#include <unistd.h>

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

typedef struct {
    unsigned char *buffer;
    size_t used;
    size_t capacity;
    int opened;
    int fd;
    int owns_fd;
    char *name;
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    MPI_File file;
    MPI_Win offset_window;
    uint64_t *offset_base;
#endif
} swbwa_output_state_t;

static swbwa_output_state_t output_state;

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
        ssize_t written = write(fd, data, length);
        if (written < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (written == 0) {
            errno = EIO;
            return -1;
        }
        data += written;
        length -= (size_t)written;
    }
    return 0;
}

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
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
    size_t position = 0;

    if (length == 0) return 0;
    increment = (uint64_t)length;
    if (increment > (uint64_t)INT64_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    if (mpi_check(MPI_Fetch_and_op(&increment, &offset, MPI_UINT64_T, 0, 0,
                                   MPI_SUM, output_state.offset_window),
                  "MPI_Fetch_and_op") != 0)
        return -1;
    if (mpi_check(MPI_Win_flush(0, output_state.offset_window),
                  "MPI_Win_flush") != 0)
        return -1;
    if (offset > (uint64_t)INT64_MAX - increment) {
        errno = EOVERFLOW;
        return -1;
    }

    while (position < length) {
        size_t remaining = length - position;
        int chunk = remaining > INT_MAX ? INT_MAX : (int)remaining;
        MPI_Status status;
        int count;

        if (mpi_check(MPI_File_write_at(output_state.file,
                                        (MPI_Offset)(offset + position),
                                        data + position, chunk, MPI_BYTE, &status),
                      "MPI_File_write_at") != 0)
            return -1;
        if (mpi_check(MPI_Get_count(&status, MPI_BYTE, &count),
                      "MPI_Get_count") != 0)
            return -1;
        if (count != chunk) {
            errno = EIO;
            return -1;
        }
        position += (size_t)chunk;
    }
    return 0;
}

static int flush_single_unordered(void)
{
    if (write_single_unordered(output_state.buffer, output_state.used) != 0)
        return -1;
    output_state.used = 0;
    return 0;
}
#endif

int swbwa_output_open(const char *path)
{
    size_t capacity = (size_t)SWBWA_OUTPUT_BUFFER_BYTES;

    if (output_state.opened) {
        errno = EALREADY;
        return -1;
    }
    if (capacity == 0 || capacity > INT_MAX) {
        errno = EINVAL;
        return -1;
    }
    memset(&output_state, 0, sizeof(output_state));
    output_state.fd = -1;
    output_state.capacity = capacity;
    output_state.buffer = malloc(capacity);
    if (output_state.buffer == NULL) return -1;

#if SWBWA_USE_MPI
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
    if (sizeof(MPI_Offset) < sizeof(int64_t)) {
        errno = EOVERFLOW;
        goto fail;
    }
    output_state.name = strdup(path);
    if (output_state.name == NULL) goto fail;
    if (mpi_check(MPI_File_open(MPI_COMM_WORLD, output_state.name,
                                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                MPI_INFO_NULL, &output_state.file),
                  "MPI_File_open") != 0)
        goto fail;
    /* MPI_File_set_size is collective; rank 0 coordinates the truncation. */
    if (mpi_check(MPI_File_set_size(output_state.file, 0),
                  "MPI_File_set_size") != 0)
        goto fail_file;
    if (mpi_check(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier") != 0)
        goto fail_file;
    if (mpi_check(MPI_Win_allocate(swbwa_mpi_rank() == 0 ? sizeof(uint64_t) : 0,
                                   sizeof(uint64_t), MPI_INFO_NULL, MPI_COMM_WORLD,
                                   &output_state.offset_base,
                                   &output_state.offset_window),
                  "MPI_Win_allocate") != 0)
        goto fail_file;
    if (mpi_check(MPI_Win_lock_all(0, output_state.offset_window),
                  "MPI_Win_lock_all") != 0)
        goto fail_window;
    if (swbwa_mpi_rank() == 0) *output_state.offset_base = 0;
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

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
fail_locked_window:
    MPI_Win_unlock_all(output_state.offset_window);
fail_window:
    MPI_Win_free(&output_state.offset_window);
fail_file:
    MPI_File_close(&output_state.file);
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
    if (length > output_state.capacity - output_state.used &&
        swbwa_output_flush() != 0)
        return -1;
    if (length > output_state.capacity) {
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
    if (!output_state.opened) {
        errno = EINVAL;
        return -1;
    }
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    return flush_single_unordered();
#else
    if (write_all(output_state.fd, output_state.buffer, output_state.used) != 0)
        return -1;
    output_state.used = 0;
    return 0;
#endif
}

int swbwa_output_close(void)
{
    int status = 0;

    if (!output_state.opened) return 0;
    if (swbwa_output_flush() != 0) status = -1;

#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_UNORDERED
    if (mpi_check(MPI_File_sync(output_state.file), "MPI_File_sync") != 0)
        status = -1;
    if (mpi_check(MPI_Win_unlock_all(output_state.offset_window),
                  "MPI_Win_unlock_all") != 0)
        status = -1;
    if (mpi_check(MPI_Win_free(&output_state.offset_window), "MPI_Win_free") != 0)
        status = -1;
    if (mpi_check(MPI_File_close(&output_state.file), "MPI_File_close") != 0)
        status = -1;
#else
    if (output_state.owns_fd && close(output_state.fd) != 0) status = -1;
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
