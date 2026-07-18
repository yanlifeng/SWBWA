#include "swbwa_config.h"
#include "swbwa_mpi.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

enum { SWBWA_MAX_BOUNDED_INPUTS = 4 };

typedef struct {
    int active;
    int fd;
    int64_t position;
    int64_t end;
} swbwa_bounded_input_t;

static swbwa_bounded_input_t bounded_inputs[SWBWA_MAX_BOUNDED_INPUTS];
static int mpi_rank;
static int mpi_size = 1;
static int mpi_initialized;

static int64_t fastq_file_size(const char *path)
{
    struct stat st;

    if (stat(path, &st) != 0) {
        fprintf(stderr, "[E::MPI input] cannot stat '%s': %s\n",
                path, strerror(errno));
        return -1;
    }
    if (!S_ISREG(st.st_mode)) {
        fprintf(stderr, "[E::MPI input] '%s' is not a regular file\n", path);
        return -1;
    }
    if (st.st_size < 0 || (uintmax_t)st.st_size > INT64_MAX) {
        fprintf(stderr, "[E::MPI input] invalid file size for '%s'\n", path);
        return -1;
    }
    return (int64_t)st.st_size;
}

static size_t fastq_line_length(const char *line, size_t length)
{
    while (length > 0 && (line[length - 1] == '\n' || line[length - 1] == '\r'))
        --length;
    return length;
}

static int fastq_record_is_valid(const char *header, ssize_t header_len,
                                 const char *sequence, ssize_t sequence_len,
                                 const char *plus, ssize_t plus_len,
                                 const char *quality, ssize_t quality_len)
{
    return header_len > 0 && header[0] == '@' &&
           sequence_len >= 0 && plus_len > 0 && plus[0] == '+' &&
           quality_len >= 0 &&
           fastq_line_length(sequence, (size_t)sequence_len) ==
           fastq_line_length(quality, (size_t)quality_len);
}

static int64_t find_fastq_boundary(FILE *file, int64_t nominal, int64_t file_size)
{
    char *header = NULL, *sequence = NULL, *plus = NULL, *quality = NULL;
    size_t header_cap = 0, sequence_cap = 0, plus_cap = 0, quality_cap = 0;
    ssize_t header_len, sequence_len, plus_len, quality_len;
    int64_t result = -1;
    int previous;

    if (nominal <= 0) return 0;
    if (nominal >= file_size) return file_size;
    if (fseeko(file, (off_t)(nominal - 1), SEEK_SET) != 0) return -1;

    previous = fgetc(file);
    if (previous == EOF) return file_size;
    if (previous != '\n' && getline(&header, &header_cap, file) < 0) {
        result = file_size;
        goto done;
    }

    for (;;) {
        off_t candidate = ftello(file);
        off_t after_header;

        if (candidate < 0) goto done;
        if ((int64_t)candidate >= file_size) {
            result = file_size;
            goto done;
        }

        header_len = getline(&header, &header_cap, file);
        if (header_len < 0) {
            result = file_size;
            goto done;
        }
        after_header = ftello(file);
        if (after_header < 0) goto done;

        if (header[0] == '@') {
            sequence_len = getline(&sequence, &sequence_cap, file);
            plus_len = getline(&plus, &plus_cap, file);
            quality_len = getline(&quality, &quality_cap, file);
            if (fastq_record_is_valid(header, header_len,
                                      sequence, sequence_len,
                                      plus, plus_len,
                                      quality, quality_len)) {
                result = (int64_t)candidate;
                goto done;
            }
        }

        if (fseeko(file, after_header, SEEK_SET) != 0) goto done;
    }

done:
    free(header);
    free(sequence);
    free(plus);
    free(quality);
    return result;
}

#if SWBWA_USE_MPI && SWBWA_ENABLE_MPI_READ_ID_SCAN
static int64_t count_fastq_records(const char *path, int64_t start, int64_t end)
{
    FILE *file = NULL;
    char *header = NULL, *sequence = NULL, *plus = NULL, *quality = NULL;
    size_t header_cap = 0, sequence_cap = 0, plus_cap = 0, quality_cap = 0;
    int64_t count = 0;

    file = fopen(path, "r");
    if (file == NULL || fseeko(file, (off_t)start, SEEK_SET) != 0) {
        fprintf(stderr,
                "[E::MPI input rank %d] cannot scan '%s' at byte %" PRId64
                ": %s\n", mpi_rank, path, start, strerror(errno));
        count = -1;
        goto done;
    }

    while (count >= 0) {
        off_t position = ftello(file);
        ssize_t header_len, sequence_len, plus_len, quality_len;

        if (position < 0) {
            count = -1;
            break;
        }
        if ((int64_t)position == end) break;
        if ((int64_t)position > end) {
            count = -1;
            break;
        }

        header_len = getline(&header, &header_cap, file);
        sequence_len = getline(&sequence, &sequence_cap, file);
        plus_len = getline(&plus, &plus_cap, file);
        quality_len = getline(&quality, &quality_cap, file);
        if (!fastq_record_is_valid(header, header_len,
                                   sequence, sequence_len,
                                   plus, plus_len,
                                   quality, quality_len)) {
            fprintf(stderr,
                    "[E::MPI input rank %d] invalid FASTQ record in '%s'"
                    " at byte %" PRId64 "\n",
                    mpi_rank, path, (int64_t)position);
            count = -1;
            break;
        }

        position = ftello(file);
        if (position < 0 || (int64_t)position > end) {
            fprintf(stderr,
                    "[E::MPI input rank %d] FASTQ range ends inside a record"
                    " in '%s' near byte %" PRId64 "\n",
                    mpi_rank, path, end);
            count = -1;
            break;
        }
        ++count;
    }

done:
    if (file != NULL) fclose(file);
    free(header);
    free(sequence);
    free(plus);
    free(quality);
    return count;
}
#endif

int swbwa_mpi_init(int *argc, char ***argv)
{
#if SWBWA_USE_MPI
    int provided;
    int rc = MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &provided);

    if (rc != MPI_SUCCESS) return -1;
    mpi_initialized = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (provided < MPI_THREAD_SERIALIZED) {
        if (mpi_rank == 0)
            fprintf(stderr, "[E::MPI] MPI_THREAD_SERIALIZED is not available\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return -1;
    }
#else
    (void)argc;
    (void)argv;
    mpi_initialized = 1;
#endif
    return 0;
}

void swbwa_mpi_finalize(void)
{
#if SWBWA_USE_MPI
    if (mpi_initialized) MPI_Finalize();
#endif
    mpi_initialized = 0;
}

int swbwa_mpi_rank(void)
{
    return mpi_rank;
}

int swbwa_mpi_size(void)
{
    return mpi_size;
}

int swbwa_mpi_is_root(void)
{
    return mpi_rank == 0;
}

void swbwa_mpi_print_rank_ordered(void (*printer)(void))
{
    if (printer == NULL) return;
#if SWBWA_USE_MPI
    {
        int rank;

        for (rank = 0; rank < mpi_size; ++rank) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (mpi_rank == rank) {
                printer();
                fflush(stderr);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    printer();
#endif
}

void swbwa_mpi_abort(const char *message)
{
    if (message != NULL)
        fprintf(stderr, "[E::MPI rank %d] %s\n", mpi_rank, message);
#if SWBWA_USE_MPI
    if (mpi_initialized) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

int swbwa_mpi_fastq_range(const char *read1_path, const char *read2_path,
                          swbwa_fastq_range_t *range)
{
#if SWBWA_USE_MPI
    int status = 0;
    int64_t *boundaries = NULL;
    int64_t file_size = 0;
#if SWBWA_ENABLE_MPI_READ_ID_SCAN
    int local_status;
    int64_t local_record_count;
    int64_t first_record = 0;
#endif

    if (mpi_rank == 0) {
        int i;
        FILE *file = NULL;

        file_size = fastq_file_size(read1_path);
        if (file_size < 0) status = -1;
        if (status == 0 && read2_path != NULL) {
            int64_t mate_size = fastq_file_size(read2_path);
            if (mate_size < 0) {
                status = -1;
            } else if (mate_size != file_size) {
                fprintf(stderr,
                        "[E::MPI input] paired FASTQ sizes differ: %" PRId64
                        " != %" PRId64 "\n", file_size, mate_size);
                status = -1;
            }
        }

        if (status == 0) {
            boundaries = calloc((size_t)mpi_size + 1, sizeof(*boundaries));
            file = fopen(read1_path, "r");
            if (boundaries == NULL || file == NULL) {
                fprintf(stderr, "[E::MPI input] cannot prepare FASTQ boundaries: %s\n",
                        strerror(errno));
                status = -1;
            }
        }

        if (status == 0) {
            boundaries[0] = 0;
            boundaries[mpi_size] = file_size;
            for (i = 1; i < mpi_size; ++i) {
                int64_t nominal = (file_size / mpi_size) * i +
                                  (file_size % mpi_size) * i / mpi_size;
                boundaries[i] = find_fastq_boundary(file, nominal, file_size);
                if (boundaries[i] < 0) {
                    fprintf(stderr,
                            "[E::MPI input] cannot locate FASTQ boundary near byte %" PRId64 "\n",
                            nominal);
                    status = -1;
                    break;
                }
                if (boundaries[i] < boundaries[i - 1])
                    boundaries[i] = boundaries[i - 1];
            }
        }
        if (file != NULL) fclose(file);
    }

    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0) {
        free(boundaries);
        return -1;
    }
    MPI_Bcast(&file_size, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        boundaries = malloc(((size_t)mpi_size + 1) * sizeof(*boundaries));
        if (boundaries == NULL) swbwa_mpi_abort("cannot allocate FASTQ boundary table");
    }
    MPI_Bcast(boundaries, mpi_size + 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    range->start = boundaries[mpi_rank];
    range->end = boundaries[mpi_rank + 1];
    range->file_size = file_size;
    free(boundaries);

#if SWBWA_ENABLE_MPI_READ_ID_SCAN
    local_record_count = count_fastq_records(read1_path, range->start, range->end);
    local_status = local_record_count < 0;
    MPI_Allreduce(&local_status, &status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (status != 0) return -1;

    MPI_Exscan(&local_record_count, &first_record, 1, MPI_INT64_T,
               MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) first_record = 0;
    range->first_record = first_record;
    range->record_count = local_record_count;
#else
    range->first_record = 0;
    range->record_count = 0;
#endif
#else
    int64_t file_size;

    (void)read2_path;
    file_size = fastq_file_size(read1_path);
    if (file_size < 0) return -1;
    range->start = 0;
    range->end = file_size;
    range->file_size = file_size;
    range->first_record = 0;
    range->record_count = 0;
#endif
    return 0;
}

int swbwa_input_register_fd(int fd, int64_t start, int64_t end)
{
    int i;

    if (start < 0 || end < start) {
        errno = EINVAL;
        return -1;
    }
#ifdef O_DIRECT
    {
        int flags = fcntl(fd, F_GETFL);

        if (flags < 0) return -1;
        /* FASTQ record boundaries do not satisfy O_DIRECT alignment. */
        if ((flags & O_DIRECT) != 0 &&
            fcntl(fd, F_SETFL, flags & ~O_DIRECT) < 0)
            return -1;
    }
#endif
    if (lseek(fd, (off_t)start, SEEK_SET) < 0) return -1;
    for (i = 0; i < SWBWA_MAX_BOUNDED_INPUTS; ++i) {
        if (!bounded_inputs[i].active) {
            bounded_inputs[i].active = 1;
            bounded_inputs[i].fd = fd;
            bounded_inputs[i].position = start;
            bounded_inputs[i].end = end;
            return 0;
        }
    }
    errno = EMFILE;
    return -1;
}

void swbwa_input_unregister_fd(int fd)
{
    int i;

    for (i = 0; i < SWBWA_MAX_BOUNDED_INPUTS; ++i) {
        if (bounded_inputs[i].active && bounded_inputs[i].fd == fd) {
            memset(&bounded_inputs[i], 0, sizeof(bounded_inputs[i]));
            return;
        }
    }
}

ssize_t swbwa_input_read(int fd, void *buffer, size_t bytes)
{
    int i;

    for (i = 0; i < SWBWA_MAX_BOUNDED_INPUTS; ++i) {
        swbwa_bounded_input_t *input = &bounded_inputs[i];
        ssize_t result;
        int64_t remaining;

        if (!input->active || input->fd != fd) continue;
        remaining = input->end - input->position;
        if (remaining <= 0) return 0;
        if ((uint64_t)remaining < bytes) bytes = (size_t)remaining;
        result = read(fd, buffer, bytes);
        if (result > 0) input->position += result;
        return result;
    }
    return read(fd, buffer, bytes);
}
