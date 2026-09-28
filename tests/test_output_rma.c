#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "swbwa_output.h"
#include "swbwa_host_workers.h"

static int rank_id, rank_count;
int swbwa_mpi_rank(void) { return rank_id; }
int swbwa_mpi_size(void) { return rank_count; }
int swbwa_mpi_is_root(void) { return rank_id == 0; }
void swbwa_mpi_print_rank_ordered(void (*printer)(void))
{
    int i;
    for (i = 0; i < rank_count; ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank_id == i) printer();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    int chunk, i, provided;
    unsigned char expected[2 * 104 * 512];
    size_t expected_bytes = 0;
    assert(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) == MPI_SUCCESS);
    assert(provided >= MPI_THREAD_FUNNELED);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
    assert(argc == 2);
    assert(swbwa_host_workers_init() == 0);
    assert(swbwa_output_open(argv[1], 1) == 0);
    /* Rank 2 deliberately has no input, but must join output teardown. */
    if (rank_id < 2) {
        for (chunk = 0; chunk < 2; ++chunk) {
            int reads = chunk ? 2 : 104;
#if TEST_CHUNK_OUTPUT
            bseq1_t seqs[104] = {{0}};
            char blobs[104][512];
            int lengths[104];
#endif
#if SWBWA_OUTPUT_RMA_ONLY && !TEST_CHUNK_OUTPUT
            assert(swbwa_output_begin_chunk(rank_id * 2 + chunk,
                                             (rank_id * 2 + chunk) * 1000,
                                             (rank_id * 2 + chunk + 1) * 1000,
                                             reads, rank_id == 1) == 0);
            assert(swbwa_output_begin_chunk(0, 0, 1, 0, 0) == -1);
            assert(swbwa_output_end_chunk() == -1);
#endif
            for (i = 0; i < reads; ++i) {
                char buffer[512];
                int length = snprintf(buffer, sizeof(buffer), "r%d-c%d-i%d\t", rank_id, chunk, i);
                /* Includes an oversized SAM blob and multi-record read output. */
                int payload = i == 0 ? 257 : 1 + i % 23;
                memset(buffer + length, 'A' + i % 20, payload);
                length += payload;
                buffer[length++] = '\n';
                if (i == 0) {
                    memcpy(buffer + length, "supplementary\n", 14);
                    length += 14;
                }
                assert(expected_bytes + (size_t)length <= sizeof(expected));
                memcpy(expected + expected_bytes, buffer, (size_t)length);
                expected_bytes += (size_t)length;
#if TEST_CHUNK_OUTPUT
                buffer[length] = '\0';
                memcpy(blobs[i], buffer, (size_t)length + 1);
                seqs[i].sam = blobs[i];
                lengths[i] = length;
#else
                assert(swbwa_output_write(buffer, length) == 0);
#endif
            }
#if TEST_CHUNK_OUTPUT
            {
                int id = rank_id * 20 + 12 - 7 * chunk;
                /* Exercise both precomputed lengths and the generic fallback. */
                assert(swbwa_output_write_chunk(id, id * 1000, (id + 1) * 1000,
                                                 seqs, chunk ? NULL : lengths, reads, rank_id == 1) == 0);
            }
#endif
#if SWBWA_OUTPUT_RMA_ONLY && !TEST_CHUNK_OUTPUT
            assert(swbwa_output_end_chunk() == 0);
            assert(swbwa_output_end_chunk() == -1);
#endif
        }
    }
    assert(swbwa_output_close() == 0);
#if !SWBWA_OUTPUT_RMA_ONLY
    /* A one-rank target-platform smoke test verifies every helper's bytes. */
    if (rank_count == 1) {
        unsigned char *actual = malloc(expected_bytes);
        FILE *file = fopen(argv[1], "rb");
        assert(actual != NULL && file != NULL);
        assert(fread(actual, 1, expected_bytes, file) == expected_bytes);
        assert(fgetc(file) == EOF);
        assert(memcmp(actual, expected, expected_bytes) == 0);
        assert(fclose(file) == 0);
        free(actual);
        printf("Full packed file: PASS (%zu bytes, all MPE slices)\n", expected_bytes);
    }
#endif
    swbwa_host_workers_destroy();
    assert(MPI_Finalize() == MPI_SUCCESS);
    return 0;
}
