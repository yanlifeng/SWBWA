#include "swbwa_config.h"
#include "swbwa_mpi.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

static int mpi_rank;
static int mpi_size = 1;
static int mpi_initialized;

#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
typedef struct {
    int64_t start;
    int64_t end;
    int64_t records;
    double stage2_seconds;
    int claimed;
    int completed;
} swbwa_chunk_debug_t;

typedef struct {
    FILE *boundary_file;
    const char *read1_path;
    int64_t *record_offsets;
    int64_t *claimed_chunk_ids;
    swbwa_chunk_debug_t *chunk_debug;
    int64_t chunk_count;
    int64_t chunk_bytes;
    int64_t file_size;
    int64_t local_chunks;
    int64_t local_records;
    int64_t local_bytes;
    int64_t next_calls;
    int64_t rma_attempts;
    int64_t remote_attempts;
    int64_t out_of_range_tickets;
    int64_t empty_chunks;
    int64_t remote_chunks;
    double setup_win_create_seconds;
    double setup_barrier_seconds;
    double next_seconds;
    double rma_seconds;
    double win_lock_seconds;
    double fetch_and_op_seconds;
    double win_unlock_seconds;
    double boundary_seconds;
    double close_barrier_seconds;
    double close_win_free_seconds;
    int debug_enabled;
    int opened;
    int next_queue;
    unsigned long long local_ticket;
    MPI_Win ticket_window;
} swbwa_fastq_scheduler_t;
#endif

typedef struct {
    char *header;
    char *sequence;
    char *plus;
    char *quality;
    size_t header_capacity;
    size_t sequence_capacity;
    size_t plus_capacity;
    size_t quality_capacity;
} swbwa_fastq_record_buffer_t;

#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
static swbwa_fastq_scheduler_t chunk_scheduler;

static double scheduler_debug_now(void)
{
    struct timespec now;

    if (!chunk_scheduler.debug_enabled) return 0.0;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec * 1.0e-9;
}

static void scheduler_debug_finish_next(double start)
{
    if (chunk_scheduler.debug_enabled)
        chunk_scheduler.next_seconds += scheduler_debug_now() - start;
}
#endif

static int64_t fastq_file_size(const char *path)
{
    struct stat st;

    if (stat(path, &st) != 0) {
        fprintf(stderr, "[E::FASTQ input] cannot stat '%s': %s\n",
                path, strerror(errno));
        return -1;
    }
    if (!S_ISREG(st.st_mode)) {
        fprintf(stderr, "[E::FASTQ input] '%s' is not a regular file\n", path);
        return -1;
    }
    if (st.st_size < 0 || (uintmax_t)st.st_size > INT64_MAX) {
        fprintf(stderr, "[E::FASTQ input] invalid file size for '%s'\n", path);
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

static void fastq_record_buffer_destroy(swbwa_fastq_record_buffer_t *buffer)
{
    free(buffer->header);
    free(buffer->sequence);
    free(buffer->plus);
    free(buffer->quality);
    memset(buffer, 0, sizeof(*buffer));
}

static int read_fastq_record(FILE *file, const char *path,
                             swbwa_fastq_record_buffer_t *buffer,
                             int64_t *start, int64_t *end,
                             int64_t *sequence_bases)
{
    off_t position;
    ssize_t header_len, sequence_len, plus_len, quality_len;

    position = ftello(file);
    if (position < 0) return -1;
    header_len = getline(&buffer->header, &buffer->header_capacity, file);
    if (header_len < 0) {
        if (feof(file)) return 0;
        fprintf(stderr, "[E::FASTQ input] cannot read '%s': %s\n",
                path, strerror(errno));
        return -1;
    }

    sequence_len = getline(&buffer->sequence, &buffer->sequence_capacity, file);
    plus_len = getline(&buffer->plus, &buffer->plus_capacity, file);
    quality_len = getline(&buffer->quality, &buffer->quality_capacity, file);
    if (!fastq_record_is_valid(buffer->header, header_len,
                               buffer->sequence, sequence_len,
                               buffer->plus, plus_len,
                               buffer->quality, quality_len)) {
        fprintf(stderr,
                "[E::FASTQ input] invalid FASTQ record in '%s' at byte %"
                PRId64 "\n",
                path, (int64_t)position);
        return -1;
    }

    *start = (int64_t)position;
    position = ftello(file);
    if (position < 0) return -1;
    *end = (int64_t)position;
    *sequence_bases = (int64_t)fastq_line_length(
        buffer->sequence, (size_t)sequence_len);
    return 1;
}

int swbwa_fastq_chunk_bytes(const char *read1_path,
                            const char *read2_path,
                            int64_t bytes_per_cg,
                            int64_t *file_size_out,
                            int64_t *chunk_bytes_out)
{
    int64_t file_size;
    int64_t chunk_bytes;

    if (bytes_per_cg <= 0 || file_size_out == NULL ||
        chunk_bytes_out == NULL) {
        errno = EINVAL;
        return -1;
    }
    file_size = fastq_file_size(read1_path);
    if (file_size < 0) return -1;
    if (read2_path != NULL) {
        int64_t mate_size = fastq_file_size(read2_path);

        if (mate_size < 0) return -1;
        if (mate_size != file_size) {
            fprintf(stderr,
                    "[E::FASTQ input] paired FASTQ sizes differ: %" PRId64
                    " != %" PRId64 "\n", file_size, mate_size);
            return -1;
        }
    }
    if (bytes_per_cg > INT64_MAX / SWBWA_CG_COUNT) {
        errno = EOVERFLOW;
        return -1;
    }

    chunk_bytes = bytes_per_cg * SWBWA_CG_COUNT;
    if (file_size == 0) chunk_bytes = 1;
    else if (chunk_bytes > file_size) chunk_bytes = file_size;
    if (chunk_bytes > SWBWA_CPE_FORMAT_BUFFER_BYTES) {
        fprintf(stderr,
                "[E::FASTQ input] chunk size %" PRId64
                " exceeds the %lld-byte CPE formatting buffer; reduce -K\n",
                chunk_bytes, (long long)SWBWA_CPE_FORMAT_BUFFER_BYTES);
        errno = E2BIG;
        return -1;
    }

    *file_size_out = file_size;
    *chunk_bytes_out = chunk_bytes;
    return 0;
}

#if SWBWA_USE_MPI
static int64_t find_fastq_boundary(FILE *file, int64_t nominal,
                                   int64_t file_size)
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
#endif

#if SWBWA_USE_MPI && SWBWA_MPI_EXACT_READ_INDEX
static int scan_fastq_record_offsets(
    FILE *file, const char *path, const int64_t *boundaries,
    int64_t boundary_count, int64_t *record_offsets)
{
    swbwa_fastq_record_buffer_t buffer = {0};
    int64_t next_boundary = 0;
    int64_t record_count = 0;
    int64_t file_end = 0;
    int status = -1;

    if (file == NULL || path == NULL || boundaries == NULL ||
        record_offsets == NULL || boundary_count <= 0) {
        errno = EINVAL;
        return -1;
    }
    for (int64_t i = 0; i < boundary_count; ++i) {
        if (boundaries[i] < 0 ||
            (i > 0 && boundaries[i] < boundaries[i - 1])) {
            errno = EINVAL;
            return -1;
        }
    }
    if (fseeko(file, 0, SEEK_SET) != 0) {
        fprintf(stderr,
                "[E::MPI input] cannot rewind '%s' for exact indexing: %s\n",
                path, strerror(errno));
        goto done;
    }

    while (next_boundary < boundary_count &&
           boundaries[next_boundary] == 0) {
        record_offsets[next_boundary++] = 0;
    }

    for (;;) {
        int64_t record_start, record_end, sequence_bases;
        int result;

        result = read_fastq_record(file, path, &buffer, &record_start,
                                   &record_end, &sequence_bases);
        if (result < 0) goto done;
        if (result == 0) break;
        (void)sequence_bases;

        if (next_boundary < boundary_count &&
            boundaries[next_boundary] < record_start) {
            fprintf(stderr,
                    "[E::MPI input] byte %" PRId64
                    " is not a FASTQ record boundary in '%s'\n",
                    boundaries[next_boundary], path);
            goto done;
        }
        while (next_boundary < boundary_count &&
               boundaries[next_boundary] == record_start) {
            record_offsets[next_boundary++] = record_count;
        }
        if (record_count == INT64_MAX) {
            errno = EOVERFLOW;
            goto done;
        }
        ++record_count;
        file_end = record_end;
    }

    while (next_boundary < boundary_count &&
           boundaries[next_boundary] == file_end) {
        record_offsets[next_boundary++] = record_count;
    }
    if (next_boundary != boundary_count) {
        fprintf(stderr,
                "[E::MPI input] exact-index boundary %" PRId64
                " lies beyond parsed FASTQ data in '%s'\n",
                boundaries[next_boundary], path);
        goto done;
    }
    status = 0;

done:
    fastq_record_buffer_destroy(&buffer);
    return status;
}

#endif

int swbwa_mpi_init(int *argc, char ***argv)
{
#if SWBWA_USE_MPI
    int provided;
    int rc = MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);

    if (rc != MPI_SUCCESS) return -1;
    mpi_initialized = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (provided < MPI_THREAD_MULTIPLE) {
        if (mpi_rank == 0)
            fprintf(stderr,
                    "[E::MPI] MPI_THREAD_MULTIPLE is required,"
                    " but MPI provided thread level %d\n",
                    provided);
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

int swbwa_mpi_barrier(void)
{
#if SWBWA_USE_MPI
    int result;

    result = MPI_Barrier(MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        errno = EIO;
        return -1;
    }
#endif
    return 0;
}

int swbwa_mpi_progress(void)
{
#if SWBWA_USE_MPI
    int flag;
    int result;
    MPI_Status status;

    result = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                        &flag, &status);
    if (result != MPI_SUCCESS) {
        errno = EIO;
        return -1;
    }
#endif
    return 0;
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
#if SWBWA_MPI_EXACT_READ_INDEX
    int64_t *record_offsets = NULL;
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
            file = fopen(read1_path, "rb");
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
#if SWBWA_MPI_EXACT_READ_INDEX
        if (status == 0) {
            record_offsets = calloc((size_t)mpi_size + 1,
                                    sizeof(*record_offsets));
            if (record_offsets == NULL ||
                scan_fastq_record_offsets(
                    file, read1_path, boundaries, mpi_size + 1,
                    record_offsets) != 0)
                status = -1;
        }
#endif
        if (file != NULL) fclose(file);
    }

    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0) {
        free(boundaries);
#if SWBWA_MPI_EXACT_READ_INDEX
        free(record_offsets);
#endif
        return -1;
    }
    MPI_Bcast(&file_size, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        boundaries = malloc(((size_t)mpi_size + 1) * sizeof(*boundaries));
        if (boundaries == NULL) swbwa_mpi_abort("cannot allocate FASTQ boundary table");
#if SWBWA_MPI_EXACT_READ_INDEX
        record_offsets = malloc(
            ((size_t)mpi_size + 1) * sizeof(*record_offsets));
        if (record_offsets == NULL)
            swbwa_mpi_abort("cannot allocate FASTQ record-offset table");
#endif
    }
    MPI_Bcast(boundaries, mpi_size + 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
#if SWBWA_MPI_EXACT_READ_INDEX
    if (MPI_Bcast(record_offsets, mpi_size + 1, MPI_INT64_T, 0,
                  MPI_COMM_WORLD) != MPI_SUCCESS) {
        free(boundaries);
        free(record_offsets);
        return -1;
    }
#endif

    range->chunk_id = -1;
    range->start = boundaries[mpi_rank];
    range->end = boundaries[mpi_rank + 1];
    range->file_size = file_size;

#if SWBWA_MPI_EXACT_READ_INDEX
    range->first_record = record_offsets[mpi_rank];
    range->record_count =
        record_offsets[mpi_rank + 1] - record_offsets[mpi_rank];
    free(record_offsets);
#else
    range->first_record = 0;
    range->record_count = 0;
#endif
    free(boundaries);
#else
    int64_t file_size;

    (void)read2_path;
    file_size = fastq_file_size(read1_path);
    if (file_size < 0) return -1;
    range->chunk_id = -1;
    range->start = 0;
    range->end = file_size;
    range->file_size = file_size;
    range->first_record = 0;
    range->record_count = 0;
#endif
    return 0;
}

#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
static int64_t chunk_nominal_offset(int64_t index)
{
    if (index <= 0) return 0;
    if (index >= chunk_scheduler.chunk_count ||
        index > chunk_scheduler.file_size / chunk_scheduler.chunk_bytes)
        return chunk_scheduler.file_size;
    return index * chunk_scheduler.chunk_bytes;
}

static int scheduler_chunk_range(int64_t index, swbwa_fastq_range_t *range)
{
    int64_t nominal_start;
    int64_t nominal_end;

    if (index < 0 || index >= chunk_scheduler.chunk_count || range == NULL) {
        errno = EINVAL;
        return -1;
    }
    nominal_start = chunk_nominal_offset(index);
    nominal_end = chunk_nominal_offset(index + 1);

    memset(range, 0, sizeof(*range));
    range->start = find_fastq_boundary(chunk_scheduler.boundary_file,
                                       nominal_start,
                                       chunk_scheduler.file_size);
    range->end = find_fastq_boundary(chunk_scheduler.boundary_file,
                                     nominal_end,
                                     chunk_scheduler.file_size);
    if (range->start < 0 || range->end < range->start) {
        fprintf(stderr,
                "[E::MPI input rank %d] cannot align FASTQ chunk %" PRId64
                " near bytes [%" PRId64 ", %" PRId64 ")\n",
                mpi_rank, index, nominal_start, nominal_end);
        return -1;
    }
    range->file_size = chunk_scheduler.file_size;
    range->chunk_id = index;
#if SWBWA_MPI_EXACT_READ_INDEX
    if (chunk_scheduler.record_offsets != NULL) {
        range->first_record = chunk_scheduler.record_offsets[index];
        range->record_count = chunk_scheduler.record_offsets[index + 1] -
                              chunk_scheduler.record_offsets[index];
    }
#endif
    return 0;
}

#if SWBWA_MPI_EXACT_READ_INDEX
static int build_exact_record_offsets(void)
{
    int64_t *boundaries = NULL;
    int value_count;
    int local_status;
    int global_status = 0;
    int root_status = 0;

    if (chunk_scheduler.chunk_count >= INT_MAX) {
        errno = EOVERFLOW;
        return -1;
    }
    value_count = (int)chunk_scheduler.chunk_count + 1;
    chunk_scheduler.record_offsets = calloc(
        (size_t)value_count,
        sizeof(*chunk_scheduler.record_offsets));
    local_status = chunk_scheduler.record_offsets == NULL;
    if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
        return -1;

    if (mpi_rank == 0 && chunk_scheduler.chunk_count > 0) {
        boundaries = malloc((size_t)value_count * sizeof(*boundaries));
        if (boundaries == NULL) {
            root_status = -1;
        } else {
            for (int index = 0; index < value_count; ++index) {
                boundaries[index] = find_fastq_boundary(
                    chunk_scheduler.boundary_file,
                    chunk_nominal_offset(index),
                    chunk_scheduler.file_size);
                if (boundaries[index] < 0 ||
                    (index > 0 &&
                     boundaries[index] < boundaries[index - 1])) {
                    root_status = -1;
                    break;
                }
            }
        }
        if (root_status == 0 &&
            scan_fastq_record_offsets(
                chunk_scheduler.boundary_file, chunk_scheduler.read1_path,
                boundaries, value_count,
                chunk_scheduler.record_offsets) != 0)
            root_status = -1;
    }
    free(boundaries);
    if (MPI_Bcast(&root_status, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
        root_status != 0)
        return -1;
    return MPI_Bcast(chunk_scheduler.record_offsets, value_count, MPI_INT64_T,
                     0, MPI_COMM_WORLD) == MPI_SUCCESS ? 0 : -1;
}
#endif

int swbwa_mpi_fastq_scheduler_open(const char *read1_path,
                                   const char *read2_path,
                                   int64_t bytes_per_cg,
                                   int debug_enabled,
                                   swbwa_fastq_range_t *assigned_range,
                                   int64_t *chunk_count)
{
    int local_status;

    if (chunk_scheduler.opened || assigned_range == NULL ||
        chunk_count == NULL || bytes_per_cg <= 0) {
        errno = EINVAL;
        return -1;
    }
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
    chunk_scheduler.read1_path = read1_path;
    chunk_scheduler.debug_enabled = debug_enabled != 0;
    local_status = swbwa_fastq_chunk_bytes(
        read1_path, read2_path, bytes_per_cg, &chunk_scheduler.file_size,
        &chunk_scheduler.chunk_bytes) != 0;

    {
        int global_status = 0;
        int64_t local_values[2];
        int64_t min_values[2];
        int64_t max_values[2];

        if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
            goto fail;
        local_values[0] = chunk_scheduler.file_size;
        local_values[1] = chunk_scheduler.chunk_bytes;
        if (MPI_Allreduce(local_values, min_values, 2, MPI_INT64_T, MPI_MIN,
                          MPI_COMM_WORLD) != MPI_SUCCESS ||
            MPI_Allreduce(local_values, max_values, 2, MPI_INT64_T, MPI_MAX,
                          MPI_COMM_WORLD) != MPI_SUCCESS)
            goto fail;
        if (min_values[0] != max_values[0] ||
            min_values[1] != max_values[1]) {
            if (mpi_rank == 0)
                fprintf(stderr,
                        "[E::MPI input] ranks disagree on FASTQ or chunk size\n");
            errno = EINVAL;
            goto fail;
        }
    }

    if (chunk_scheduler.file_size != 0) {
        chunk_scheduler.chunk_count =
            chunk_scheduler.file_size / chunk_scheduler.chunk_bytes +
            (chunk_scheduler.file_size % chunk_scheduler.chunk_bytes != 0);
        chunk_scheduler.boundary_file = fopen(read1_path, "rb");
        local_status = chunk_scheduler.boundary_file == NULL;
    } else {
        local_status = 0;
    }
    if (chunk_scheduler.debug_enabled && chunk_scheduler.chunk_count > 0) {
        if ((uint64_t)chunk_scheduler.chunk_count >
            SIZE_MAX / sizeof(*chunk_scheduler.claimed_chunk_ids) ||
            (uint64_t)chunk_scheduler.chunk_count >
            SIZE_MAX / sizeof(*chunk_scheduler.chunk_debug)) {
            errno = EOVERFLOW;
            local_status = 1;
        } else {
            chunk_scheduler.claimed_chunk_ids = malloc(
                (size_t)chunk_scheduler.chunk_count *
                sizeof(*chunk_scheduler.claimed_chunk_ids));
            chunk_scheduler.chunk_debug = calloc(
                (size_t)chunk_scheduler.chunk_count,
                sizeof(*chunk_scheduler.chunk_debug));
            if (chunk_scheduler.claimed_chunk_ids == NULL ||
                chunk_scheduler.chunk_debug == NULL)
                local_status = 1;
        }
    }

    {
        int global_status = 0;

        if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
            goto fail;
    }

#if SWBWA_MPI_EXACT_READ_INDEX
    if (build_exact_record_offsets() != 0) goto fail;
#endif

    chunk_scheduler.local_ticket = 0;
    chunk_scheduler.next_queue = mpi_rank;
    {
        double start = scheduler_debug_now();

        if (MPI_Win_create(&chunk_scheduler.local_ticket,
                           sizeof(chunk_scheduler.local_ticket),
                           sizeof(chunk_scheduler.local_ticket), MPI_INFO_NULL,
                           MPI_COMM_WORLD,
                           &chunk_scheduler.ticket_window) != MPI_SUCCESS)
            goto fail;
        if (chunk_scheduler.debug_enabled)
            chunk_scheduler.setup_win_create_seconds +=
                scheduler_debug_now() - start;
    }
    {
        double start = scheduler_debug_now();

        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
            goto fail_window;
        if (chunk_scheduler.debug_enabled)
            chunk_scheduler.setup_barrier_seconds +=
                scheduler_debug_now() - start;
    }

    memset(assigned_range, 0, sizeof(*assigned_range));
    assigned_range->chunk_id = -1;
    assigned_range->file_size = chunk_scheduler.file_size;
    assigned_range->start = 0;
    assigned_range->end = chunk_scheduler.file_size;
#if SWBWA_MPI_EXACT_READ_INDEX
    if (chunk_scheduler.record_offsets != NULL)
        assigned_range->record_count =
            chunk_scheduler.record_offsets[chunk_scheduler.chunk_count];
#endif
    *chunk_count = chunk_scheduler.chunk_count;
    chunk_scheduler.opened = 1;
    return 0;

fail_window:
    MPI_Win_free(&chunk_scheduler.ticket_window);
fail:
    if (chunk_scheduler.boundary_file != NULL)
        fclose(chunk_scheduler.boundary_file);
    free(chunk_scheduler.record_offsets);
    free(chunk_scheduler.claimed_chunk_ids);
    free(chunk_scheduler.chunk_debug);
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
    return -1;
}

int64_t swbwa_mpi_fastq_scheduler_chunk_bytes(void)
{
    return chunk_scheduler.chunk_bytes;
}

int swbwa_mpi_fastq_scheduler_next(swbwa_fastq_range_t *range)
{
    double next_start;

    if (!chunk_scheduler.opened || range == NULL) {
        errno = EINVAL;
        return -1;
    }
    next_start = scheduler_debug_now();
    if (chunk_scheduler.debug_enabled) ++chunk_scheduler.next_calls;

    for (;;) {
        /* Queue q owns chunk IDs q, q + mpi_size, ... . */
        int start_queue = chunk_scheduler.next_queue;
        int claimed_valid_chunk = 0;
        int attempt;

        for (attempt = 0; attempt < mpi_size; ++attempt) {
            const int queue = (start_queue + attempt) % mpi_size;
            const unsigned long long increment = 1;
            unsigned long long ticket = 0;
            int64_t index;
            int result;
            int unlock_result;
            int lock_succeeded;
            double rma_start = scheduler_debug_now();
            double operation_start;

            if (chunk_scheduler.debug_enabled) {
                ++chunk_scheduler.rma_attempts;
                if (queue != mpi_rank) ++chunk_scheduler.remote_attempts;
            }
            operation_start = scheduler_debug_now();
            result = MPI_Win_lock(MPI_LOCK_SHARED, queue, 0,
                                  chunk_scheduler.ticket_window);
            lock_succeeded = result == MPI_SUCCESS;
            if (chunk_scheduler.debug_enabled) {
                double operation_end = scheduler_debug_now();

                chunk_scheduler.win_lock_seconds +=
                    operation_end - operation_start;
                operation_start = operation_end;
            }
            if (result == MPI_SUCCESS) {
                result = MPI_Fetch_and_op(
                    &increment, &ticket, MPI_UNSIGNED_LONG_LONG, queue, 0,
                    MPI_SUM, chunk_scheduler.ticket_window);
                if (chunk_scheduler.debug_enabled)
                    chunk_scheduler.fetch_and_op_seconds +=
                        scheduler_debug_now() - operation_start;
            }
            if (lock_succeeded) {
                operation_start = scheduler_debug_now();
                unlock_result = MPI_Win_unlock(
                    queue, chunk_scheduler.ticket_window);
                if (chunk_scheduler.debug_enabled)
                    chunk_scheduler.win_unlock_seconds +=
                        scheduler_debug_now() - operation_start;
                if (result == MPI_SUCCESS) result = unlock_result;
            }
            if (chunk_scheduler.debug_enabled)
                chunk_scheduler.rma_seconds +=
                    scheduler_debug_now() - rma_start;

            if (result != MPI_SUCCESS) {
                errno = EIO;
                scheduler_debug_finish_next(next_start);
                return -1;
            }
            if (ticket >
                (unsigned long long)(INT64_MAX - queue) /
                    (unsigned long long)mpi_size) {
                errno = EOVERFLOW;
                scheduler_debug_finish_next(next_start);
                return -1;
            }

            index = queue + (int64_t)ticket * mpi_size;
            if (index >= chunk_scheduler.chunk_count) {
                if (chunk_scheduler.debug_enabled)
                    ++chunk_scheduler.out_of_range_tickets;
                continue;
            }

            claimed_valid_chunk = 1;
            chunk_scheduler.next_queue = queue;
            {
                double boundary_start = scheduler_debug_now();
                int range_status = scheduler_chunk_range(index, range);

                if (chunk_scheduler.debug_enabled)
                    chunk_scheduler.boundary_seconds +=
                        scheduler_debug_now() - boundary_start;
                if (range_status != 0) {
                    scheduler_debug_finish_next(next_start);
                    return -1;
                }
            }
            if (range->start == range->end) {
                if (chunk_scheduler.debug_enabled)
                    ++chunk_scheduler.empty_chunks;
                break;
            }

            if (chunk_scheduler.debug_enabled) {
                swbwa_chunk_debug_t *debug_chunk =
                    &chunk_scheduler.chunk_debug[index];

                chunk_scheduler.claimed_chunk_ids[
                    chunk_scheduler.local_chunks] = index;
                debug_chunk->start = range->start;
                debug_chunk->end = range->end;
                debug_chunk->claimed = 1;
                if (queue != mpi_rank) ++chunk_scheduler.remote_chunks;
            }
            ++chunk_scheduler.local_chunks;
            chunk_scheduler.local_bytes += range->end - range->start;
            scheduler_debug_finish_next(next_start);
            return 1;
        }
        if (!claimed_valid_chunk) {
            scheduler_debug_finish_next(next_start);
            return 0;
        }
    }
}

void swbwa_mpi_fastq_scheduler_record_stage2(int64_t chunk_id,
                                             int64_t records,
                                             double seconds)
{
    if (records > 0) chunk_scheduler.local_records += records;
    if (chunk_scheduler.debug_enabled &&
        chunk_id >= 0 && chunk_id < chunk_scheduler.chunk_count) {
        swbwa_chunk_debug_t *debug_chunk =
            &chunk_scheduler.chunk_debug[chunk_id];

        if (!debug_chunk->claimed) return;
        debug_chunk->records = records;
        debug_chunk->stage2_seconds = seconds;
        debug_chunk->completed = 1;
    }
}

void swbwa_mpi_fastq_scheduler_stats(int64_t *chunks, int64_t *records,
                                     int64_t *bytes)
{
    if (chunks != NULL) *chunks = chunk_scheduler.local_chunks;
    if (records != NULL) *records = chunk_scheduler.local_records;
    if (bytes != NULL) *bytes = chunk_scheduler.local_bytes;
}

static void print_scheduler_debug_time(const char *label, double seconds,
                                       int64_t calls, int indent)
{
    fprintf(stderr, "%*s%-42s %10.6f s", indent, "", label, seconds);
    if (calls > 0)
        fprintf(stderr, "  (%" PRId64 " calls, %9.3f us/call)",
                calls, seconds * 1.0e6 / calls);
    fputc('\n', stderr);
}

static void print_scheduler_debug_report_body(void)
{
    double local_next_seconds =
        chunk_scheduler.next_seconds -
        chunk_scheduler.rma_seconds -
        chunk_scheduler.boundary_seconds;
    double other_rma_seconds =
        chunk_scheduler.rma_seconds -
        chunk_scheduler.win_lock_seconds -
        chunk_scheduler.fetch_and_op_seconds -
        chunk_scheduler.win_unlock_seconds;
    double stage2_seconds = 0.0;
    double max_stage2_seconds = 0.0;
    int64_t completed_chunks = 0;
    int64_t max_stage2_chunk = -1;
    int64_t order;

    if (local_next_seconds < 0.0) local_next_seconds = 0.0;
    if (other_rma_seconds < 0.0) other_rma_seconds = 0.0;

    for (order = 0; order < chunk_scheduler.local_chunks; ++order) {
        int64_t chunk_id = chunk_scheduler.claimed_chunk_ids[order];
        const swbwa_chunk_debug_t *debug_chunk =
            &chunk_scheduler.chunk_debug[chunk_id];

        if (!debug_chunk->completed) continue;
        ++completed_chunks;
        stage2_seconds += debug_chunk->stage2_seconds;
        if (max_stage2_chunk < 0 ||
            debug_chunk->stage2_seconds > max_stage2_seconds) {
            max_stage2_seconds = debug_chunk->stage2_seconds;
            max_stage2_chunk = chunk_id;
        }
    }

    fprintf(stderr,
            "\n"
            "================ MPI Dynamic Scheduler Debug =================\n"
            "  MPI rank: %06d / %06d\n"
            "\n"
            "  Work summary\n"
            "    claimed chunks                               %12" PRId64 "\n"
            "    claimed from this rank's queue               %12" PRId64 "\n"
            "    claimed from another rank's queue            %12" PRId64 "\n"
            "    completed FASTQ records                      %12" PRId64 "\n"
            "    claimed input bytes                          %12" PRId64 "\n"
            "\n"
            "  Scheduler counters\n"
            "    scheduler_next calls                         %12" PRId64 "\n"
            "    RMA ticket attempts                          %12" PRId64 "\n"
            "    attempts against another rank's queue        %12" PRId64 "\n"
            "    tickets beyond the final chunk               %12" PRId64 "\n"
            "    empty aligned chunks                         %12" PRId64 "\n",
            mpi_rank, mpi_size,
            chunk_scheduler.local_chunks,
            chunk_scheduler.local_chunks - chunk_scheduler.remote_chunks,
            chunk_scheduler.remote_chunks,
            chunk_scheduler.local_records,
            chunk_scheduler.local_bytes,
            chunk_scheduler.next_calls,
            chunk_scheduler.rma_attempts,
            chunk_scheduler.remote_attempts,
            chunk_scheduler.out_of_range_tickets,
            chunk_scheduler.empty_chunks);

    fprintf(stderr, "\n  Scheduler timing\n");
    print_scheduler_debug_time(
        "MPI_Win_create (setup)",
        chunk_scheduler.setup_win_create_seconds, 1, 4);
    print_scheduler_debug_time(
        "MPI_Barrier after setup",
        chunk_scheduler.setup_barrier_seconds, 1, 4);
    print_scheduler_debug_time(
        "scheduler_next total",
        chunk_scheduler.next_seconds, chunk_scheduler.next_calls, 4);
    print_scheduler_debug_time(
        "RMA ticket path total",
        chunk_scheduler.rma_seconds, chunk_scheduler.rma_attempts, 6);
    print_scheduler_debug_time(
        "MPI_Win_lock",
        chunk_scheduler.win_lock_seconds,
        chunk_scheduler.rma_attempts, 8);
    print_scheduler_debug_time(
        "MPI_Fetch_and_op",
        chunk_scheduler.fetch_and_op_seconds,
        chunk_scheduler.rma_attempts, 8);
    print_scheduler_debug_time(
        "MPI_Win_unlock",
        chunk_scheduler.win_unlock_seconds,
        chunk_scheduler.rma_attempts, 8);
    print_scheduler_debug_time(
        "RMA call overhead (derived)",
        other_rma_seconds, chunk_scheduler.rma_attempts, 8);
    print_scheduler_debug_time(
        "align nominal offsets to FASTQ records",
        chunk_scheduler.boundary_seconds,
        chunk_scheduler.local_chunks + chunk_scheduler.empty_chunks, 6);
    print_scheduler_debug_time(
        "scheduler loop/bookkeeping (derived)",
        local_next_seconds, chunk_scheduler.next_calls, 6);

    fprintf(stderr, "\n  Scheduler teardown timing\n");
    print_scheduler_debug_time(
        "MPI_Barrier (wait for other ranks)",
        chunk_scheduler.close_barrier_seconds, 1, 4);
    print_scheduler_debug_time(
        "MPI_Win_free",
        chunk_scheduler.close_win_free_seconds, 1, 4);

    fprintf(stderr,
            "\n"
            "  Stage 2 by chunk\n"
            "    completed chunks                             %12" PRId64 "\n"
            "    accumulated stage 2                          %12.6f s\n"
            "    average stage 2 per chunk                    %12.6f s\n"
            "    slowest chunk                                %12" PRId64 "\n"
            "    slowest chunk stage 2                        %12.6f s\n"
            "\n"
            "    order   chunk  queue        byte range"
            "                       bytes      records   stage2(s)\n",
            completed_chunks,
            stage2_seconds,
            completed_chunks > 0 ? stage2_seconds / completed_chunks : 0.0,
            max_stage2_chunk,
            max_stage2_seconds);

    for (order = 0; order < chunk_scheduler.local_chunks; ++order) {
        int64_t chunk_id = chunk_scheduler.claimed_chunk_ids[order];
        const swbwa_chunk_debug_t *debug_chunk =
            &chunk_scheduler.chunk_debug[chunk_id];

        fprintf(stderr,
                "    %5" PRId64 " %7" PRId64 " %6d"
                "  [%12" PRId64 ", %12" PRId64 ")"
                " %12" PRId64 " %12" PRId64 " %11.6f%s\n",
                order, chunk_id, (int)(chunk_id % mpi_size),
                debug_chunk->start, debug_chunk->end,
                debug_chunk->end - debug_chunk->start,
                debug_chunk->records,
                debug_chunk->stage2_seconds,
                debug_chunk->completed ? "" : "  incomplete");
    }
    fprintf(stderr,
            "==============================================================\n"
            "\n");
}

static void scheduler_debug_report(void)
{
    if (!chunk_scheduler.opened || !chunk_scheduler.debug_enabled) return;
    swbwa_mpi_print_rank_ordered(print_scheduler_debug_report_body);
}

void swbwa_mpi_fastq_scheduler_close(void)
{
    double start;

    if (!chunk_scheduler.opened) return;

    start = scheduler_debug_now();
    MPI_Barrier(MPI_COMM_WORLD);
    if (chunk_scheduler.debug_enabled)
        chunk_scheduler.close_barrier_seconds +=
            scheduler_debug_now() - start;

    start = scheduler_debug_now();
    MPI_Win_free(&chunk_scheduler.ticket_window);
    if (chunk_scheduler.debug_enabled)
        chunk_scheduler.close_win_free_seconds +=
            scheduler_debug_now() - start;

    scheduler_debug_report();
    if (chunk_scheduler.boundary_file != NULL)
        fclose(chunk_scheduler.boundary_file);
    free(chunk_scheduler.record_offsets);
    free(chunk_scheduler.claimed_chunk_ids);
    free(chunk_scheduler.chunk_debug);
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
}
#endif
