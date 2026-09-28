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

int main(int argc, char **argv)
{
    int rank, records, pass;
    assert(swbwa_mpi_init(&argc, &argv) == 0);
    assert(argc == 5);
    rank = swbwa_mpi_rank();
    records = atoi(argv[4]);
    assert(records >= 0 && records <= 100);
    assert(swbwa_host_workers_init() == 0);
    for (pass = 0; pass < 2; ++pass) {
        int64_t positions[101], count;
        swbwa_fastq_range_t initial, range;
        int result, i;
        if (rank == 0) {
            FILE *fp = fopen(argv[1], "wb");
            assert(fp);
            for (i = 0; i < records; ++i) {
                positions[i] = ftello(fp);
                fprintf(fp, "@read%d\nACGTACGT\n+\nIIIIIIII\n", i);
            }
            positions[records] = ftello(fp);
            assert(fclose(fp) == 0);
        }
        MPI_Bcast(positions, records + 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        assert(swbwa_mpi_fastq_scheduler_open(argv[1], NULL, atoll(argv[3]),
                                             1, &initial, &count) == 0);
        assert(strcmp(swbwa_mpi_fastq_scheduler_ticket_mode(), "global") == 0);
        assert(swbwa_output_open(argv[2], 1) == 0);
        assert(swbwa_output_write("unidentified", 12) == -1 && errno == EINVAL);
        assert(swbwa_mpi_progress_thread_start(0) == 0);
        while ((result = swbwa_mpi_fastq_scheduler_next(&range)) > 0) {
            bseq1_t seqs[100] = {{0}};
            char sam[100][400];
            int lengths[100], n = 0;
            for (i = 0; i < records; ++i) {
                if (positions[i] < range.start || positions[i] >= range.end) continue;
                lengths[n] = snprintf(sam[n], sizeof(sam[n]),
                    "read%03d\t0\tchr1\t%d\nread%03d\t2048\tchr2\t%d\n", i, i+1, i, i+2);
                seqs[n].sam = sam[n];
                ++n;
            }
            assert(n > 0);
            if (range.chunk_id == 0) {
                struct timespec delay = {0, 20000000};
                nanosleep(&delay, NULL);
            }
            assert(swbwa_output_write_chunk(range.chunk_id, range.start, range.end,
                                            seqs, lengths, n, 0) == 0);
            swbwa_mpi_fastq_scheduler_record_stage2(range.chunk_id, n, 0);
        }
        assert(result == 0);
        assert(swbwa_output_close() == 0);
        MPI_Barrier(MPI_COMM_WORLD);
#if !SWBWA_OUTPUT_RMA_ONLY
        if (rank == 0) {
            FILE *fp = fopen(argv[2], "rb");
            char expected[400], actual[400];
            assert(fp);
            for (i = 0; i < records; ++i) {
                int n = snprintf(expected, sizeof(expected),
                    "read%03d\t0\tchr1\t%d\nread%03d\t2048\tchr2\t%d\n", i, i+1, i, i+2);
                assert(fread(actual, 1, (size_t)n, fp) == (size_t)n);
                assert(memcmp(actual, expected, (size_t)n) == 0);
            }
            assert(fgetc(fp) == EOF);
            assert(fclose(fp) == 0);
            puts("ordered output file bytes PASS");
        }
#endif
        swbwa_mpi_fastq_scheduler_close();
        assert(swbwa_mpi_progress_thread_stop() == 0);
    }
    swbwa_host_workers_destroy();
    swbwa_mpi_finalize();
    if (rank == 0) puts("ordered output PASS");
    return 0;
}
