#include "swbwa_mpi.h"
#include "swbwa_output.h"
#include "swbwa_host_workers.h"

#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

enum { MAX_READS = 137, FINISH_TAG = 73 };

static int chunk_reads(int pass, int chunk)
{
    return chunk % 3 == 0 ? MAX_READS - pass * 30 : 5 + chunk % 5 - pass;
}

static char *sam_blob(int pass, int chunk, int read, int zero_bytes)
{
    size_t capacity = read == 100 ? 72000 : 1024, used = 0;
    char *sam = calloc(capacity, 1);
    int record;
    assert(sam);
    if (zero_bytes || chunk % 4 == 1 || read % 7 == 6) return sam;
    for (record = 0; record < 1 + (chunk + read) % 3; ++record) {
        int n = snprintf(sam + used, capacity - used,
            "p%d-c%03d-r%03d\t%d\tchr%d\t%d\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII",
            pass, chunk, read, record ? 2048 : 0, record + 1, chunk * 1000 + read + 1);
        assert(n > 0 && (size_t)n < capacity - used);
        used += (size_t)n;
        /* A blob beyond the default 100 reads also exceeds 64 KiB. */
        if (read == 100 && record == 0) {
            size_t padding = 70000 - (size_t)pass * 123;
            memcpy(sam + used, "\tZZ:Z:", 6);
            used += 6;
            assert(used + padding + 2 < capacity);
            memset(sam + used, 'A' + chunk % 20, padding);
            used += padding;
            sam[used++] = 'a' + pass;
        }
        sam[used++] = '\n';
        sam[used] = '\0';
    }
    return sam;
}

static void write_chunk(int pass, int chunk, int zero_bytes)
{
    bseq1_t seqs[MAX_READS] = {{0}};
    int lengths[MAX_READS], i, reads = chunk_reads(pass, chunk);
    for (i = 0; i < reads; ++i) {
        seqs[i].sam = sam_blob(pass, chunk, i, zero_bytes);
        lengths[i] = (int)strlen(seqs[i].sam);
    }
    assert(swbwa_output_write_chunk(chunk, chunk * 10000LL,
                                    (chunk + 1) * 10000LL, seqs,
                                    pass ? NULL : lengths, reads, 0) == 0);
    for (i = 0; i < reads; ++i) free(seqs[i].sam);
}

static void check_file(const char *path, int pass, int chunks, int zero_bytes)
{
#if SWBWA_OUTPUT_RMA_ONLY
    (void)pass;
    (void)chunks;
    (void)zero_bytes;
    errno = 0;
    assert(access(path, F_OK) == -1 && errno == ENOENT);
#else
    FILE *file = fopen(path, "rb");
    int chunk, read;
    assert(file);
    for (chunk = 0; chunk < chunks; ++chunk) {
        for (read = 0; read < chunk_reads(pass, chunk); ++read) {
            char *expected = sam_blob(pass, chunk, read, zero_bytes);
            size_t length = strlen(expected);
            char *actual = malloc(length + 1);
            assert(actual);
            assert(fread(actual, 1, length, file) == length);
            assert(memcmp(actual, expected, length) == 0);
            free(actual);
            free(expected);
        }
    }
    assert(fgetc(file) == EOF && !ferror(file));
    assert(fclose(file) == 0);
#endif
}

int main(int argc, char **argv)
{
    int rank, ranks, chunks, pass, zero_bytes, decreasing, missing, missing_one;
    int wrong_low, wrong_high, inconsistent, expect_error;
    assert(swbwa_mpi_init(&argc, &argv) == 0);
    assert(argc == 4);
    rank = swbwa_mpi_rank();
    ranks = swbwa_mpi_size();
    chunks = atoi(argv[2]);
    assert(chunks >= 0 && chunks <= 100);
    zero_bytes = strcmp(argv[3], "zero-bytes") == 0;
    decreasing = strcmp(argv[3], "decreasing") == 0;
    missing = strcmp(argv[3], "missing-finish") == 0;
    missing_one = strcmp(argv[3], "missing-one-finish") == 0;
    wrong_low = strcmp(argv[3], "wrong-low") == 0;
    wrong_high = strcmp(argv[3], "wrong-high") == 0;
    inconsistent = strcmp(argv[3], "inconsistent-count") == 0;
    expect_error = missing || missing_one || wrong_low || wrong_high || inconsistent;
    assert(strcmp(argv[3], "roundrobin") == 0 || zero_bytes || decreasing || expect_error);
    assert(!wrong_low || chunks > 0);
    assert(!decreasing || chunks > ranks);
    assert(swbwa_host_workers_init() == 0);

    for (pass = 0; pass < 2; ++pass) {
        int chunk, status, failed, failures, ack = 1;
        int local_finish_check = !expect_error && ranks > 1 && chunks < ranks;
        int64_t final_chunks = chunks;
        if (rank == 0) {
            fprintf(stderr, "[online test] pass=%d begin\n", pass);
            fflush(stderr);
        }
        assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
        /* No FASTQ scheduler or chunk-count publication precedes this open. */
        assert(swbwa_output_open_online(argv[1], 1) == 0);
        errno = 0;
        assert(swbwa_output_write("unidentified", 12) == -1 && errno == EINVAL);
        assert(swbwa_mpi_progress_thread_start(0) == 0);

        /* An idle rank must finish locally while rank 0 has not written yet. */
        if (local_finish_check && rank == 0) {
            assert(MPI_Recv(&ack, 1, MPI_INT, ranks - 1, FINISH_TAG,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE) == MPI_SUCCESS);
            assert(ack == 1);
            puts("online local finish PASS");
        }
        if (rank == 0 && chunks) {
            struct timespec delay = {0, 50000000};
            while (nanosleep(&delay, &delay) != 0) assert(errno == EINTR);
        }
        for (chunk = rank; chunk < chunks; chunk += ranks)
            write_chunk(pass, chunk, zero_bytes);
        if (decreasing && rank == 0) {
            int previous = (chunks - 1) / ranks * ranks - ranks;
            errno = 0;
            assert(swbwa_output_write_chunk(previous, 0, 0, NULL, NULL, 0, 0) == -1);
            assert(errno == EINVAL);
            puts("online decreasing ID rejected PASS");
        }
        assert(swbwa_output_flush() == 0);
        if (wrong_low) --final_chunks;
        if (wrong_high || (inconsistent && rank == ranks - 1)) ++final_chunks;
        /* EOF is known and this synchronous writer has completely drained. */
        if (!missing && !(missing_one && rank == ranks - 1))
            assert(swbwa_output_ordered_finish(final_chunks) == 0);
        if (local_finish_check && rank == ranks - 1)
            assert(MPI_Send(&ack, 1, MPI_INT, 0, FINISH_TAG, MPI_COMM_WORLD) == MPI_SUCCESS);

        status = swbwa_output_close();
        assert(status == (expect_error ? -1 : 0));
        failed = status != 0;
        assert(MPI_Allreduce(&failed, &failures, 1, MPI_INT, MPI_SUM,
                              MPI_COMM_WORLD) == MPI_SUCCESS);
        assert(failures == (expect_error ? ranks : 0));
        assert(swbwa_mpi_progress_thread_stop() == 0);
        if (rank == 0) {
            check_file(argv[1], pass, chunks, zero_bytes);
            printf("online %s pass=%d PASS\n", expect_error ? "close error" : "output bytes", pass);
            fprintf(stderr, "[online test] pass=%d end\n", pass);
            fflush(stderr);
        }
        assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
    }
    swbwa_host_workers_destroy();
    swbwa_mpi_finalize();
    if (rank == 0) puts("online ordered output PASS");
    return 0;
}
