#include "swbwa_config.h"
#include "swbwa_mpi.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "malloc_wrap.h"

#if SWBWA_USE_MPI
#include <mpi.h>
#endif

enum {
    SWBWA_MAX_BOUNDED_INPUTS = 4,
    SWBWA_FASTQ_ESTIMATE_RECORDS = 1024
};

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
#if SWBWA_USE_MPI
static pthread_mutex_t mpi_call_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#if SWBWA_ENABLE_MPI_FASTQ_SCHEDULER
typedef struct {
    FILE *boundary_file;
    const char *read1_path;
    int64_t *record_offsets;
    int64_t chunk_count;
    int64_t chunk_bytes;
    int64_t file_size;
    int64_t local_chunks;
    int64_t local_records;
    int64_t local_bytes;
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

#if SWBWA_ENABLE_MPI_FASTQ_SCHEDULER
static swbwa_fastq_scheduler_t chunk_scheduler;
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

int swbwa_fastq_estimate_chunk_bytes(const char *read1_path,
                                     const char *read2_path,
                                     int64_t target_bases,
                                     int64_t *file_size_out,
                                     int64_t *chunk_bytes_out)
{
    FILE *read1 = NULL;
    FILE *read2 = NULL;
    swbwa_fastq_record_buffer_t buffer1 = {0};
    swbwa_fastq_record_buffer_t buffer2 = {0};
    int64_t file_size;
    int64_t sampled_bytes = 0;
    int64_t sampled_bases = 0;
    int64_t max_record_bytes = 0;
    int sampled_records = 0;
    int status = -1;

    if (target_bases <= 0 || file_size_out == NULL || chunk_bytes_out == NULL) {
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

    if (file_size == 0) {
        *file_size_out = 0;
        *chunk_bytes_out = 1;
        return 0;
    }

    read1 = fopen(read1_path, "rb");
    if (read2_path != NULL) read2 = fopen(read2_path, "rb");
    if (read1 == NULL || (read2_path != NULL && read2 == NULL)) {
        fprintf(stderr, "[E::FASTQ input] cannot open input: %s\n",
                strerror(errno));
        goto done;
    }

    while (sampled_records < SWBWA_FASTQ_ESTIMATE_RECORDS) {
        int result1, result2 = 0;
        int64_t start1, end1, bases1;
        int64_t start2 = 0, end2 = 0, bases2 = 0;
        int64_t record_bytes;

        result1 = read_fastq_record(read1, read1_path, &buffer1,
                                    &start1, &end1, &bases1);
        if (read2 != NULL)
            result2 = read_fastq_record(read2, read2_path, &buffer2,
                                        &start2, &end2, &bases2);
        if (result1 < 0 || result2 < 0) goto done;
        if (result1 == 0 || (read2 != NULL && result2 == 0)) {
            if (result1 != 0 || (read2 != NULL && result2 != 0)) {
                fprintf(stderr,
                        "[E::FASTQ input] paired FASTQ record counts differ"
                        " in the initial sample\n");
                goto done;
            }
            break;
        }
        if (read2 != NULL && (start1 != start2 || end1 != end2)) {
            fprintf(stderr,
                    "[E::FASTQ input] paired FASTQ byte layouts differ in the"
                    " initial sample; identical offsets are required\n");
            goto done;
        }
        if (bases1 > INT64_MAX - bases2 ||
            sampled_bases > INT64_MAX - bases1 - bases2 ||
            sampled_bytes > INT64_MAX - (end1 - start1)) {
            errno = EOVERFLOW;
            goto done;
        }
        record_bytes = end1 - start1;
        sampled_bytes += record_bytes;
        sampled_bases += bases1 + bases2;
        if (record_bytes > max_record_bytes)
            max_record_bytes = record_bytes;
        ++sampled_records;
    }

    if (sampled_records == 0 || sampled_bases == 0 || sampled_bytes == 0) {
        fprintf(stderr,
                "[E::FASTQ input] cannot estimate chunk size from '%s'\n",
                read1_path);
        errno = EINVAL;
        goto done;
    }

    {
        long double estimate = (long double)target_bases * sampled_bytes /
                               sampled_bases;
        int64_t chunk_bytes;

        if (estimate > INT64_MAX) {
            errno = EOVERFLOW;
            goto done;
        }
        chunk_bytes = (int64_t)estimate;
        if ((long double)chunk_bytes < estimate) ++chunk_bytes;
        if (chunk_bytes < max_record_bytes)
            chunk_bytes = max_record_bytes;
        if (chunk_bytes <= 0) chunk_bytes = 1;
        if (chunk_bytes > file_size) chunk_bytes = file_size;
        *file_size_out = file_size;
        *chunk_bytes_out = chunk_bytes;
    }
    status = 0;

done:
    if (read1 != NULL) fclose(read1);
    if (read2 != NULL) fclose(read2);
    fastq_record_buffer_destroy(&buffer1);
    fastq_record_buffer_destroy(&buffer2);
    return status;
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

#if SWBWA_USE_MPI && SWBWA_ENABLE_MPI_EXACT_READ_INDEX
static int64_t count_fastq_records_in_range(
    FILE *file, const char *path, int64_t start, int64_t end,
    swbwa_fastq_record_buffer_t *buffer)
{
    int64_t count = 0;

    if (fseeko(file, (off_t)start, SEEK_SET) != 0) {
        fprintf(stderr,
                "[E::MPI input rank %d] cannot scan '%s' at byte %" PRId64
                ": %s\n", mpi_rank, path, start, strerror(errno));
        return -1;
    }

    while (1) {
        int result;
        int64_t record_start, record_end, sequence_bases;
        off_t position = ftello(file);

        if (position < 0) return -1;
        if ((int64_t)position == end) break;
        if ((int64_t)position > end) return -1;

        result = read_fastq_record(file, path, buffer, &record_start,
                                   &record_end, &sequence_bases);
        if (result <= 0 || record_start != (int64_t)position ||
            record_end > end) {
            fprintf(stderr,
                    "[E::MPI input rank %d] invalid FASTQ record in '%s'"
                    " at byte %" PRId64 "\n",
                    mpi_rank, path, (int64_t)position);
            return -1;
        }
        ++count;
    }
    return count;
}

static int64_t count_fastq_records(const char *path, int64_t start, int64_t end)
{
    FILE *file = fopen(path, "rb");
    swbwa_fastq_record_buffer_t buffer = {0};
    int64_t count;

    if (file == NULL) {
        fprintf(stderr, "[E::MPI input rank %d] cannot open '%s': %s\n",
                mpi_rank, path, strerror(errno));
        return -1;
    }
    count = count_fastq_records_in_range(file, path, start, end, &buffer);
    fclose(file);
    fastq_record_buffer_destroy(&buffer);
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

void swbwa_mpi_call_lock(void)
{
#if SWBWA_USE_MPI
    pthread_mutex_lock(&mpi_call_mutex);
#endif
}

void swbwa_mpi_call_unlock(void)
{
#if SWBWA_USE_MPI
    pthread_mutex_unlock(&mpi_call_mutex);
#endif
}

int swbwa_mpi_barrier(void)
{
#if SWBWA_USE_MPI
    int result;

    swbwa_mpi_call_lock();
    result = MPI_Barrier(MPI_COMM_WORLD);
    swbwa_mpi_call_unlock();
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
#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
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

#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
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

#if SWBWA_ENABLE_MPI_FASTQ_SCHEDULER
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
#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
    if (chunk_scheduler.record_offsets != NULL) {
        range->first_record = chunk_scheduler.record_offsets[index];
        range->record_count = chunk_scheduler.record_offsets[index + 1] -
                              chunk_scheduler.record_offsets[index];
    }
#endif
    return 0;
}

#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
static int allreduce_int64_sum(int64_t *values, int64_t count)
{
    int64_t offset = 0;

    while (offset < count) {
        int block = count - offset > INT_MAX
                  ? INT_MAX : (int)(count - offset);

        if (MPI_Allreduce(MPI_IN_PLACE, values + offset, block, MPI_INT64_T,
                          MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
            return -1;
        offset += block;
    }
    return 0;
}

static int build_exact_record_offsets(void)
{
    swbwa_fastq_record_buffer_t buffer = {0};
    int64_t index;
    int64_t first_chunk;
    int64_t local_chunk_count;
    int local_status = 0;
    int global_status = 0;

    if (chunk_scheduler.chunk_count >
        (int64_t)(SIZE_MAX / sizeof(*chunk_scheduler.record_offsets)) - 1) {
        errno = EOVERFLOW;
        return -1;
    }
    chunk_scheduler.record_offsets = calloc(
        (size_t)chunk_scheduler.chunk_count + 1,
        sizeof(*chunk_scheduler.record_offsets));
    local_status = chunk_scheduler.record_offsets == NULL;
    if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
        return -1;
    local_status = 0;
    global_status = 0;

    local_chunk_count = chunk_scheduler.chunk_count / mpi_size;
    first_chunk = local_chunk_count * mpi_rank;
    if (mpi_rank < chunk_scheduler.chunk_count % mpi_size) {
        ++local_chunk_count;
        first_chunk += mpi_rank;
    } else {
        first_chunk += chunk_scheduler.chunk_count % mpi_size;
    }

    for (index = first_chunk; index < first_chunk + local_chunk_count;
         ++index) {
        swbwa_fastq_range_t range;
        int64_t count;

        if (scheduler_chunk_range(index, &range) != 0) {
            local_status = 1;
            break;
        }
        count = count_fastq_records_in_range(
            chunk_scheduler.boundary_file, chunk_scheduler.read1_path,
            range.start, range.end, &buffer);
        if (count < 0) {
            local_status = 1;
            break;
        }
        chunk_scheduler.record_offsets[index + 1] = count;
    }
    fastq_record_buffer_destroy(&buffer);

    if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
        return -1;
    if (allreduce_int64_sum(chunk_scheduler.record_offsets,
                            chunk_scheduler.chunk_count + 1) != 0)
        return -1;

    for (index = 1; index <= chunk_scheduler.chunk_count; ++index) {
        if (chunk_scheduler.record_offsets[index] >
            INT64_MAX - chunk_scheduler.record_offsets[index - 1]) {
            errno = EOVERFLOW;
            return -1;
        }
        chunk_scheduler.record_offsets[index] +=
            chunk_scheduler.record_offsets[index - 1];
    }
    return 0;
}
#endif

int swbwa_mpi_fastq_scheduler_open(const char *read1_path,
                                   const char *read2_path,
                                   int64_t target_bases,
                                   swbwa_fastq_range_t *assigned_range,
                                   int64_t *chunk_count)
{
    int local_status;

    if (chunk_scheduler.opened || assigned_range == NULL ||
        chunk_count == NULL || target_bases <= 0) {
        errno = EINVAL;
        return -1;
    }
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
    chunk_scheduler.read1_path = read1_path;
    local_status = swbwa_fastq_estimate_chunk_bytes(
        read1_path, read2_path, target_bases, &chunk_scheduler.file_size,
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
                        "[E::MPI input] ranks disagree on FASTQ size or chunk estimate\n");
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

    {
        int global_status = 0;

        if (MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD) != MPI_SUCCESS || global_status != 0)
            goto fail;
    }

#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
    if (build_exact_record_offsets() != 0) goto fail;
#endif

    chunk_scheduler.local_ticket = 0;
    chunk_scheduler.next_queue = mpi_rank;
    if (MPI_Win_create(&chunk_scheduler.local_ticket,
                       sizeof(chunk_scheduler.local_ticket),
                       sizeof(chunk_scheduler.local_ticket), MPI_INFO_NULL,
                       MPI_COMM_WORLD,
                       &chunk_scheduler.ticket_window) != MPI_SUCCESS)
        goto fail;
    if (MPI_Win_lock_all(0, chunk_scheduler.ticket_window) != MPI_SUCCESS)
        goto fail_window;
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        goto fail_locked_window;

    memset(assigned_range, 0, sizeof(*assigned_range));
    assigned_range->file_size = chunk_scheduler.file_size;
    assigned_range->start = 0;
    assigned_range->end = chunk_scheduler.file_size;
#if SWBWA_ENABLE_MPI_EXACT_READ_INDEX
    if (chunk_scheduler.record_offsets != NULL)
        assigned_range->record_count =
            chunk_scheduler.record_offsets[chunk_scheduler.chunk_count];
#endif
    *chunk_count = chunk_scheduler.chunk_count;
    chunk_scheduler.opened = 1;
    return 0;

fail_locked_window:
    MPI_Win_unlock_all(chunk_scheduler.ticket_window);
fail_window:
    MPI_Win_free(&chunk_scheduler.ticket_window);
fail:
    if (chunk_scheduler.boundary_file != NULL)
        fclose(chunk_scheduler.boundary_file);
    free(chunk_scheduler.record_offsets);
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
    return -1;
}

int64_t swbwa_mpi_fastq_scheduler_chunk_bytes(void)
{
    return chunk_scheduler.chunk_bytes;
}

int swbwa_mpi_fastq_scheduler_next(swbwa_fastq_range_t *range)
{
    if (!chunk_scheduler.opened || range == NULL) {
        errno = EINVAL;
        return -1;
    }
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

            swbwa_mpi_call_lock();
            result = MPI_Fetch_and_op(
                &increment, &ticket, MPI_UNSIGNED_LONG_LONG, queue, 0,
                MPI_SUM, chunk_scheduler.ticket_window);
            if (result == MPI_SUCCESS)
                result = MPI_Win_flush(queue, chunk_scheduler.ticket_window);
            swbwa_mpi_call_unlock();
            if (result != MPI_SUCCESS) {
                errno = EIO;
                return -1;
            }
            if (ticket >
                (unsigned long long)(INT64_MAX - queue) /
                    (unsigned long long)mpi_size) {
                errno = EOVERFLOW;
                return -1;
            }

            index = queue + (int64_t)ticket * mpi_size;
            if (index >= chunk_scheduler.chunk_count) continue;

            claimed_valid_chunk = 1;
            chunk_scheduler.next_queue = queue;
            if (scheduler_chunk_range(index, range) != 0) return -1;
            if (range->start == range->end) break;

            ++chunk_scheduler.local_chunks;
            chunk_scheduler.local_bytes += range->end - range->start;
            return 1;
        }
        if (!claimed_valid_chunk) return 0;
    }
}

void swbwa_mpi_fastq_scheduler_add_records(int64_t records)
{
    if (records > 0) chunk_scheduler.local_records += records;
}

void swbwa_mpi_fastq_scheduler_stats(int64_t *chunks, int64_t *records,
                                     int64_t *bytes)
{
    if (chunks != NULL) *chunks = chunk_scheduler.local_chunks;
    if (records != NULL) *records = chunk_scheduler.local_records;
    if (bytes != NULL) *bytes = chunk_scheduler.local_bytes;
}

void swbwa_mpi_fastq_scheduler_close(void)
{
    if (!chunk_scheduler.opened) return;
    MPI_Win_flush_all(chunk_scheduler.ticket_window);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_unlock_all(chunk_scheduler.ticket_window);
    MPI_Win_free(&chunk_scheduler.ticket_window);
    if (chunk_scheduler.boundary_file != NULL)
        fclose(chunk_scheduler.boundary_file);
    free(chunk_scheduler.record_offsets);
    memset(&chunk_scheduler, 0, sizeof(chunk_scheduler));
}
#endif

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

int swbwa_input_set_range(int fd, int64_t start, int64_t end)
{
    int i;

    if (start < 0 || end < start) {
        errno = EINVAL;
        return -1;
    }
    for (i = 0; i < SWBWA_MAX_BOUNDED_INPUTS; ++i) {
        swbwa_bounded_input_t *input = &bounded_inputs[i];

        if (!input->active || input->fd != fd) continue;
        if (lseek(fd, (off_t)start, SEEK_SET) < 0) return -1;
        input->position = start;
        input->end = end;
        return 0;
    }
    errno = ENOENT;
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
