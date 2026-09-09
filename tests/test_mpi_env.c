#include <mpi.h>

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

enum { TEST_RECORD_BYTES = 64 };

static const char *thread_level_name(int level)
{
    switch (level) {
    case MPI_THREAD_SINGLE:
        return "MPI_THREAD_SINGLE";
    case MPI_THREAD_FUNNELED:
        return "MPI_THREAD_FUNNELED";
    case MPI_THREAD_SERIALIZED:
        return "MPI_THREAD_SERIALIZED";
    case MPI_THREAD_MULTIPLE:
        return "MPI_THREAD_MULTIPLE";
    default:
        return "unknown";
    }
}

static void mpi_fatal(int result, const char *operation, int rank)
{
    char message[MPI_MAX_ERROR_STRING];
    int length = 0;

    if (result == MPI_SUCCESS) return;
    MPI_Error_string(result, message, &length);
    fprintf(stderr, "[rank %d] %s failed: %.*s\n",
            rank, operation, length, message);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    exit(EXIT_FAILURE);
}

static void make_record(char record[TEST_RECORD_BYTES], int rank,
                        uint64_t offset)
{
    char text[TEST_RECORD_BYTES];
    int length;

    memset(record, ' ', TEST_RECORD_BYTES);
    length = snprintf(text, sizeof(text),
                      "rank=%06d offset=%020" PRIu64, rank, offset);
    if (length < 0) length = 0;
    if (length > TEST_RECORD_BYTES - 1) length = TEST_RECORD_BYTES - 1;
    memcpy(record, text, (size_t)length);
    record[TEST_RECORD_BYTES - 1] = '\n';
}

int main(int argc, char **argv)
{
    const char *output_path = argc > 1 ? argv[1] : "test_mpi_env.dat";
    const uint64_t increment = TEST_RECORD_BYTES;
    uint64_t claimed_offset = UINT64_MAX;
    uint64_t expected_size;
    uint64_t final_offset = 0;
    uint64_t *counter = NULL;
    MPI_Win counter_window;
    MPI_File output_file;
    MPI_Status write_status;
    char record[TEST_RECORD_BYTES];
    int requested = MPI_THREAD_SERIALIZED;
    int provided = MPI_THREAD_SINGLE;
    int mpi_major = 0, mpi_minor = 0;
    int rank = 0, size = 1;
    int write_count = 0;
    int local_error = 0, global_error = 0;
    int i;

    if (MPI_Init_thread(&argc, &argv, requested, &provided) != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init_thread failed\n");
        return EXIT_FAILURE;
    }
    mpi_fatal(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "MPI_Comm_rank", rank);
    mpi_fatal(MPI_Comm_size(MPI_COMM_WORLD, &size), "MPI_Comm_size", rank);
    mpi_fatal(MPI_Get_version(&mpi_major, &mpi_minor), "MPI_Get_version", rank);

    if (rank == 0) {
        printf("============================================================\n");
        printf("MPI environment test\n");
        printf("  MPI standard version : %d.%d\n", mpi_major, mpi_minor);
        printf("  ranks                : %d\n", size);
        printf("  requested thread     : %s (%d)\n",
               thread_level_name(requested), requested);
        printf("  provided thread      : %s (%d)\n",
               thread_level_name(provided), provided);
        if (provided < MPI_THREAD_SERIALIZED) {
            printf("  thread result        : insufficient for MPI calls from the\n");
            printf("                         current SWBWA writer thread\n");
        } else {
            printf("  thread result        : sufficient\n");
        }
        printf("============================================================\n");
        fflush(stdout);
    }

    /* Keep all MPI calls below on the thread that called MPI_Init_thread. */
    mpi_fatal(MPI_Win_allocate(rank == 0 ? (MPI_Aint)sizeof(uint64_t) : 0,
                               (int)sizeof(uint64_t), MPI_INFO_NULL,
                               MPI_COMM_WORLD, &counter, &counter_window),
              "MPI_Win_allocate", rank);
    mpi_fatal(MPI_Win_lock_all(0, counter_window), "MPI_Win_lock_all", rank);
    if (rank == 0) *counter = 0;
    mpi_fatal(MPI_Win_sync(counter_window), "MPI_Win_sync(init)", rank);
    mpi_fatal(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(init)", rank);

    mpi_fatal(MPI_Fetch_and_op(&increment, &claimed_offset, MPI_UINT64_T,
                               0, 0, MPI_SUM, counter_window),
              "MPI_Fetch_and_op", rank);
    mpi_fatal(MPI_Win_flush(0, counter_window), "MPI_Win_flush", rank);
    mpi_fatal(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(RMA)", rank);
    mpi_fatal(MPI_Win_sync(counter_window), "MPI_Win_sync(result)", rank);

    expected_size = (uint64_t)size * TEST_RECORD_BYTES;
    if (rank == 0) final_offset = *counter;
    mpi_fatal(MPI_Bcast(&final_offset, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD),
              "MPI_Bcast(final offset)", rank);

    if (claimed_offset % TEST_RECORD_BYTES != 0 ||
        claimed_offset >= expected_size)
        local_error = 1;
    if (final_offset != expected_size) local_error = 1;

    mpi_fatal(MPI_Win_unlock_all(counter_window), "MPI_Win_unlock_all", rank);

    for (i = 0; i < size; ++i) {
        mpi_fatal(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(print)", rank);
        if (rank == i) {
            printf("[rank %d] claimed RMA byte range [%-6" PRIu64 ", %-6" PRIu64 ")\n",
                   rank, claimed_offset, claimed_offset + increment);
            fflush(stdout);
        }
    }
    mpi_fatal(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(after print)", rank);

    mpi_fatal(MPI_File_open(MPI_COMM_WORLD, (char *)output_path,
                            MPI_MODE_CREATE | MPI_MODE_WRONLY,
                            MPI_INFO_NULL, &output_file),
              "MPI_File_open", rank);
    /* MPI_File_set_size is collective, so every rank must call it. */
    mpi_fatal(MPI_File_set_size(output_file, 0), "MPI_File_set_size", rank);
    mpi_fatal(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(file truncate)", rank);

    make_record(record, rank, claimed_offset);
    mpi_fatal(MPI_File_write_at(output_file, (MPI_Offset)claimed_offset,
                                record, TEST_RECORD_BYTES, MPI_BYTE,
                                &write_status),
              "MPI_File_write_at", rank);
    mpi_fatal(MPI_Get_count(&write_status, MPI_BYTE, &write_count),
              "MPI_Get_count", rank);
    if (write_count != TEST_RECORD_BYTES) local_error = 1;

    mpi_fatal(MPI_File_sync(output_file), "MPI_File_sync", rank);
    mpi_fatal(MPI_File_close(&output_file), "MPI_File_close", rank);
    mpi_fatal(MPI_Win_free(&counter_window), "MPI_Win_free", rank);

    if (rank == 0) {
        struct stat st;

        if (stat(output_path, &st) != 0) {
            fprintf(stderr, "stat('%s') failed: %s\n",
                    output_path, strerror(errno));
            local_error = 1;
        } else if ((uint64_t)st.st_size != expected_size) {
            fprintf(stderr,
                    "unexpected file size: actual=%" PRIu64
                    " expected=%" PRIu64 "\n",
                    (uint64_t)st.st_size, expected_size);
            local_error = 1;
        }
    }

    mpi_fatal(MPI_Allreduce(&local_error, &global_error, 1, MPI_INT,
                            MPI_MAX, MPI_COMM_WORLD),
              "MPI_Allreduce", rank);
    if (rank == 0) {
        printf("------------------------------------------------------------\n");
        printf("RMA Fetch-and-op test : %s\n",
               final_offset == expected_size ? "PASS" : "FAIL");
        printf("MPI-IO write test     : %s\n",
               global_error == 0 ? "PASS" : "FAIL");
        printf("Output file           : %s (%" PRIu64 " bytes expected)\n",
               output_path, expected_size);
        printf("Overall               : %s\n",
               global_error == 0 ? "PASS" : "FAIL");
        printf("------------------------------------------------------------\n");
        fflush(stdout);
    }

    if (MPI_Finalize() != MPI_SUCCESS) {
        fprintf(stderr, "[rank %d] MPI_Finalize failed\n", rank);
        return EXIT_FAILURE;
    }
    return global_error == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
