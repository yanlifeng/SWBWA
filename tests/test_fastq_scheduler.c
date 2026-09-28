#include "swbwa_mpi.h"
#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv)
{
    int rank, pass;
    assert(swbwa_mpi_init(&argc, &argv) == 0);
    rank = swbwa_mpi_rank();
    assert(argc == 4 || argc == 5);
    if (argc == 5) {
        int64_t chunks;
        swbwa_fastq_range_t range;
        if (strcmp(argv[4], "invalid") == 0)
            setenv("SWBWA_MPI_TICKET_MODE", "invalid", 1);
        else {
            assert(strcmp(argv[4], "mixed") == 0);
            setenv("SWBWA_MPI_TICKET_MODE", rank == 1 ? "global" : "distributed", 1);
        }
        assert(swbwa_mpi_fastq_scheduler_open(argv[1], NULL,
                atoll(argv[2]), 1, &range, &chunks) == -1);
        assert(errno == EINVAL);
        swbwa_mpi_finalize();
        if (rank == 0) puts("FASTQ scheduler policy rejection PASS");
        return 0;
    }
    for (pass = 0; pass < 2; ++pass) {
        int records = atoi(argv[3]), i, result;
        int64_t positions[101], chunks, *local, *all;
        swbwa_fastq_range_t initial, range;
        assert(records >= 0 && records <= 100);
        if (rank == 0) {
            FILE *fp = fopen(argv[1], "wb");
            assert(fp != NULL);
            for (i = 0; i < records; ++i) {
                int j, length = 31 + i % 5;
                positions[i] = ftello(fp);
                fprintf(fp, "@read%d\n", i);
                for (j = 0; j < length; ++j) fputc("ACGT"[j % 4], fp);
                fputs("\n+\n", fp);
                for (j = 0; j < length; ++j) fputc('I', fp);
                fputc('\n', fp);
            }
            positions[records] = ftello(fp);
            assert(fclose(fp) == 0);
        }
        MPI_Bcast(positions, records + 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        assert(swbwa_mpi_fastq_scheduler_open(argv[1], argv[1],
                atoll(argv[2]), 1, &initial, &chunks) == 0);
        local = calloc((size_t)(chunks + 1) * 4, sizeof(*local));
        all = calloc((size_t)(chunks + 1) * 4, sizeof(*all));
        assert(local != NULL && all != NULL);
        while ((result = swbwa_mpi_fastq_scheduler_next(&range)) > 0) {
            int first = 0, end = 0;
            assert(range.chunk_id >= 0 && range.chunk_id < chunks);
            for (i = 0; i <= records; ++i) {
                if (positions[i] == range.start) first = i;
                if (positions[i] == range.end) end = i;
            }
            assert(end > first && positions[first] == range.start);
#if SWBWA_MPI_EXACT_READ_INDEX
            assert(range.first_record == first && range.record_count == end - first);
#endif
            assert(local[range.chunk_id * 4] == 0);
            local[range.chunk_id * 4] = 1;
            local[range.chunk_id * 4 + 1] = range.start;
            local[range.chunk_id * 4 + 2] = range.end;
            local[range.chunk_id * 4 + 3] = end - first;
            swbwa_mpi_fastq_scheduler_record_stage2(range.chunk_id, end - first, 0);
            if (rank == 1) {
                struct timespec delay = {0, 100000};
                nanosleep(&delay, NULL);
            }
        }
        assert(result == 0);
        assert(swbwa_mpi_fastq_scheduler_next(&range) == 0);
        MPI_Allreduce(local, all, (int)(chunks + 1) * 4, MPI_INT64_T,
                      MPI_SUM, MPI_COMM_WORLD);
        {
            int64_t end = 0, n = 0;
            for (i = 0; i < chunks; ++i) {
                assert(all[i * 4] <= 1);
                if (!all[i * 4]) continue; /* Empty aligned chunks. */
                assert(all[i * 4 + 1] == end);
                end = all[i * 4 + 2];
                n += all[i * 4 + 3];
            }
            assert(end == positions[records] && n == records);
        }
        swbwa_mpi_fastq_scheduler_close();
        free(all);
        free(local);
    }
    swbwa_mpi_finalize();
    if (rank == 0) puts("FASTQ scheduler coverage PASS");
    return 0;
}
