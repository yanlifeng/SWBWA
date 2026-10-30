#include "swbwa_output.h"
#include "swbwa_mpi.h"
#include "swbwa_host_workers.h"
#include "codec_job.h"
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

static pthread_t producer;
static int prepare_api, expect_prepared, empty_records, fail_codec;
static void *encoder_workspace;
static int64_t known_chunks;

int64_t __wrap_swbwa_mpi_fastq_scheduler_chunk_count(void) { return known_chunks; }
const char *__wrap_swbwa_mpi_fastq_scheduler_ticket_mode(void) { return "global"; }

void __real_swbwa_codec_launch(swbwa_codec_batch *batch);
void __wrap_swbwa_codec_launch(swbwa_codec_batch *batch)
{
    if (expect_prepared) assert(pthread_equal(pthread_self(), producer));
    if (encoder_workspace) assert(batch->workspaces == encoder_workspace);
    else encoder_workspace = batch->workspaces;
    __real_swbwa_codec_launch(batch);
    if (fail_codec) { batch->jobs[0].status = -1; fail_codec = 0; }
}

typedef struct {
    int chunk, reads;
    bseq1_t seqs[31];
    int lengths[31];
    void *prepared;
} output_batch;

static void free_sam(output_batch *batch)
{
    int read;
    for (read = 0; read < batch->reads; ++read) {
        free(batch->seqs[read].sam);
        batch->seqs[read].sam = NULL;
    }
}

static output_batch *prepare(int chunk, int pass)
{
    output_batch *batch = calloc(1, sizeof(*batch));
    bseq1_t *seqs;
    int *lengths, read;
    assert(batch);
    batch->chunk = chunk;
    batch->reads = empty_records ? 0 : 31;
    seqs = batch->seqs;
    lengths = batch->lengths;
    for (read = 0; read < batch->reads; ++read) {
        size_t bytes = (read * 73 + chunk * 41) % 4000;
        char prefix[128];
        int n = snprintf(prefix, sizeof(prefix), "chunk%03d-read%03d-pass%d\t", chunk, read, pass);
        size_t i;
        if (read == 5 && chunk % 2 == 0) bytes += 1048576;
        if (read % 7 == 0) bytes = 0;
        lengths[read] = bytes ? (int)(bytes + (size_t)n + 1) : 0;
        seqs[read].sam = calloc((size_t)lengths[read] + 1, 1);
        assert(seqs[read].sam);
        if (bytes) {
            memcpy(seqs[read].sam, prefix, (size_t)n);
            for (i = 0; i < bytes; ++i) {
                if (chunk == 2 && read == 5) {
                    uint32_t mixed = (uint32_t)i * UINT32_C(2654435761);
                    mixed ^= mixed >> 13;
                    mixed *= UINT32_C(0x85ebca6b);
                    mixed ^= mixed >> 16;
                    seqs[read].sam[n + i] = (char)(33 + mixed % 90);
                } else seqs[read].sam[n + i] = "ACGTN0123456789"[(i * 7 + read + chunk) % 15];
            }
            seqs[read].sam[lengths[read] - 1] = '\n';
        }
    }
    if (prepare_api) {
        assert(!swbwa_output_prepare_chunk(seqs, chunk % 2 ? NULL : lengths,
                                           batch->reads, &batch->prepared));
        assert((batch->prepared != NULL) == expect_prepared);
        assert(!swbwa_output_publish_prepared_chunk(chunk, batch->reads, batch->prepared));
#if SWBWA_USE_MPI && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (batch->prepared && getenv("SWBWA_GZIP_PREFIX_AHEAD") &&
            !strcmp(getenv("SWBWA_GZIP_PREFIX_AHEAD"), "1")) {
            assert(swbwa_output_publish_prepared_chunk(chunk, batch->reads, batch->prepared) == -1 && errno == EINVAL);
            assert(swbwa_output_release_prepared_chunk(batch->prepared) == -1 && errno == EBUSY);
        }
#endif
#if !SWBWA_OUTPUT_RMA_ONLY
        if (batch->prepared) free_sam(batch);
#endif
    }
    return batch;
}

static void consume(output_batch *batch)
{
    int chunk = batch->chunk;
    const bseq1_t *seqs = batch->seqs;
#if !SWBWA_OUTPUT_RMA_ONLY
    if (batch->prepared) seqs = NULL;
#endif
    if (prepare_api)
        assert(!swbwa_output_write_prepared_chunk(chunk, chunk * 1000, (chunk + 1) * 1000,
                    seqs, chunk % 2 ? NULL : batch->lengths, batch->reads, 0, batch->prepared));
    else
        assert(!swbwa_output_write_chunk(chunk, chunk * 1000, (chunk + 1) * 1000,
                    seqs, chunk % 2 ? NULL : batch->lengths, batch->reads, 0));
    free_sam(batch);
    free(batch);
}

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t changed;
    output_batch **batches;
    int produced, consumed, done, initial;
} queue;

static void *writer(void *opaque)
{
    queue *q = opaque;
    assert(!pthread_mutex_lock(&q->mutex));
    while (q->produced < q->initial && !q->done)
        assert(!pthread_cond_wait(&q->changed, &q->mutex));
    for (;;) {
        output_batch *batch;
        while (q->produced == q->consumed && !q->done)
            assert(!pthread_cond_wait(&q->changed, &q->mutex));
        if (q->produced == q->consumed) break;
        batch = q->batches[q->consumed];
        assert(!pthread_mutex_unlock(&q->mutex));
        usleep(10000);
        consume(batch);
        assert(!pthread_mutex_lock(&q->mutex));
        ++q->consumed;
        assert(!pthread_cond_broadcast(&q->changed));
    }
    assert(!pthread_mutex_unlock(&q->mutex));
    return NULL;
}

static void emit_parallel(int rank, int ranks, int chunks, int pass)
{
    pthread_t thread;
    int chunk;
    queue q = {.mutex = PTHREAD_MUTEX_INITIALIZER, .changed = PTHREAD_COND_INITIALIZER};
    q.batches = calloc((size_t)chunks + 1, sizeof(*q.batches));
    assert(q.batches);
    q.initial = SWBWA_PIPELINE_BUFFER_COUNT;
    assert(!pthread_create(&thread, NULL, writer, &q));
    for (chunk = rank; chunk < chunks; chunk += ranks) {
        output_batch *batch;
        assert(!pthread_mutex_lock(&q.mutex));
        while (q.produced - q.consumed == SWBWA_PIPELINE_BUFFER_COUNT)
            assert(!pthread_cond_wait(&q.changed, &q.mutex));
        assert(!pthread_mutex_unlock(&q.mutex));
        batch = prepare(chunk, pass);
        assert(!pthread_mutex_lock(&q.mutex));
        q.batches[q.produced++] = batch;
        assert(!pthread_cond_broadcast(&q.changed));
        assert(!pthread_mutex_unlock(&q.mutex));
    }
    assert(!pthread_mutex_lock(&q.mutex));
    q.done = 1;
    assert(!pthread_cond_broadcast(&q.changed));
    assert(!pthread_mutex_unlock(&q.mutex));
    assert(!pthread_join(thread, NULL));
    assert(q.produced == q.consumed);
    assert(!pthread_cond_destroy(&q.changed));
    assert(!pthread_mutex_destroy(&q.mutex));
    free(q.batches);
}

static void check_pool_errors(void)
{
    void *held[SWBWA_PIPELINE_BUFFER_COUNT], *extra = (void *)1;
    bseq1_t invalid = {0};
    int i, negative = -1;
    assert(expect_prepared);
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &held[i]) && held[i]);
    assert(swbwa_output_prepare_chunk(NULL, NULL, 0, &extra) == -1 && errno == EAGAIN && !extra);
    assert(swbwa_output_close() == -1 && errno == EBUSY);
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        assert(!swbwa_output_release_prepared_chunk(held[i]));
    assert(swbwa_output_release_prepared_chunk(held[0]) == -1 && errno == EINVAL);
    assert(!swbwa_output_release_prepared_chunk(NULL));
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT * 2; ++i) {
        assert(swbwa_output_prepare_chunk(&invalid, NULL, 1, &extra) == -1 && errno == EINVAL && !extra);
        invalid.sam = "invalid\n";
        assert(swbwa_output_prepare_chunk(&invalid, &negative, 1, &extra) == -1 && errno == EINVAL && !extra);
        invalid.sam = NULL;
        fail_codec = 1;
        assert(swbwa_output_prepare_chunk(NULL, NULL, 0, &extra) == -1 && errno == EIO && !extra);
        assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &extra) && extra);
        assert(swbwa_output_write_prepared_chunk(-1, 0, 0, NULL, NULL, 0, 0, extra) == -1 && errno == EINVAL);
        assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &extra) && extra);
        assert(swbwa_output_write_prepared_chunk(0, 0, 0, &invalid, NULL, 1, 0, extra) == -1 && errno == EINVAL);
    }
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &held[i]) && held[i]);
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        assert(!swbwa_output_release_prepared_chunk(held[i]));
}

int main(int argc, char **argv)
{
    int rank, ranks, pass, chunks;
    const char *prepare_setting = getenv("SWBWA_GZIP_PREPARE_OUTPUT");
    assert(swbwa_mpi_init(&argc, &argv) == 0);
    assert(argc == 3);
    chunks = atoi(argv[2]);
    known_chunks = chunks;
    rank = swbwa_mpi_rank();
    ranks = swbwa_mpi_size();
    producer = pthread_self();
    prepare_api = getenv("SWBWA_TEST_PREPARE_API") != NULL;
    expect_prepared = prepare_api && (!prepare_setting || !strcmp(prepare_setting, "1")) &&
                      strlen(argv[1]) >= 3 && !strcmp(argv[1] + strlen(argv[1]) - 3, ".gz");
    empty_records = getenv("SWBWA_TEST_EMPTY_RECORDS") != NULL;
    if (!getenv("SWBWA_TEST_NO_HOST_POOL")) assert(swbwa_host_workers_init() == 0);
    for (pass = 0; pass < 2; ++pass) {
        int chunk;
        encoder_workspace = NULL;
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (getenv("SWBWA_TEST_KNOWN_PREFIX"))
            assert(swbwa_output_open(argv[1], 1) == 0);
        else assert(swbwa_output_open_online(argv[1], 1) == 0);
#else
        assert(swbwa_output_open(argv[1], 1) == 0);
#endif
#if SWBWA_USE_MPI
        assert(swbwa_mpi_progress_thread_start(0) == 0);
#endif
        if (getenv("SWBWA_TEST_POOL_ERRORS")) check_pool_errors();
        if (getenv("SWBWA_TEST_WRITE_FAILURE")) {
            void *held[SWBWA_PIPELINE_BUFFER_COUNT];
            assert(!SWBWA_USE_MPI && expect_prepared);
            for (chunk = 0; chunk < SWBWA_PIPELINE_BUFFER_COUNT * 2; ++chunk) {
                void *handle;
                assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &handle) && handle);
                assert(swbwa_output_write_prepared_chunk(chunk, 0, 0, NULL, NULL, 0, 0, handle) == -1 && errno == ENOSPC);
            }
            for (chunk = 0; chunk < SWBWA_PIPELINE_BUFFER_COUNT; ++chunk)
                assert(!swbwa_output_prepare_chunk(NULL, NULL, 0, &held[chunk]));
            for (chunk = 0; chunk < SWBWA_PIPELINE_BUFFER_COUNT; ++chunk)
                assert(!swbwa_output_release_prepared_chunk(held[chunk]));
            assert(swbwa_output_close() == -1);
            continue;
        }
        if (getenv("SWBWA_TEST_SLOW_WRITER")) emit_parallel(rank, ranks, chunks, pass);
        else for (chunk = rank; chunk < chunks; chunk += ranks) consume(prepare(chunk, pass));
#if SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_SINGLE_ORDERED
        if (!getenv("SWBWA_TEST_KNOWN_PREFIX"))
            assert(swbwa_output_ordered_finish(chunks) == 0);
#endif
        assert(swbwa_output_close() == 0);
#if SWBWA_USE_MPI
        assert(swbwa_mpi_progress_thread_stop() == 0);
#endif
    }
    swbwa_host_workers_destroy();
    swbwa_mpi_finalize();
    if (!rank) puts("gzip output driver PASS");
    return 0;
}
