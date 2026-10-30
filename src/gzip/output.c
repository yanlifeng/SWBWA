#include "swbwa_gzip_output.h"
#include "swbwa_config.h"
#include "swbwa_host_workers.h"
#include "codec_job.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct swbwa_gzip_output_batch {
    struct iovec *vectors;
    unsigned char *slab;
    size_t slab_capacity;
    unsigned vector_capacity, count;
    uint64_t raw_bytes, gzip_bytes;
    int reads, in_use;
};

struct swbwa_gzip_output {
    swbwa_codec_span *spans;
    size_t *prefix;
    size_t record_capacity;
    swbwa_codec_job *jobs;
    unsigned char *workspace_storage, *workspaces;
    unsigned job_capacity, workers;
    pthread_mutex_t encoder_mutex, slots_mutex;
    swbwa_gzip_output_batch slots[SWBWA_PIPELINE_BUFFER_COUNT], *legacy;
    size_t block_bytes;
    int level;
    const bseq1_t *seqs;
    const int *lengths;
    int reads, errors[SWBWA_HOST_MPE_THREADS];
    size_t sums[SWBWA_HOST_MPE_THREADS], offsets[SWBWA_HOST_MPE_THREADS];
    uint64_t batches, members, input_bytes, output_bytes, ldm_members;
    unsigned max_active_workers;
    double prepare_seconds, cpe_seconds, wait_seconds;
};

static void mutex_check(int error)
{
    if (!error) return;
    fprintf(stderr, "[E::gzip output] synchronization failed: %s\n", strerror(error));
    abort();
}

static int encoder_lock(swbwa_gzip_output *writer)
{
    int error = pthread_mutex_trylock(&writer->encoder_mutex);
    if (!error) return 0;
    errno = error;
    return -1;
}

static swbwa_gzip_output_batch *acquire_batch(swbwa_gzip_output *writer)
{
    unsigned i;
    swbwa_gzip_output_batch *batch = NULL;
    mutex_check(pthread_mutex_lock(&writer->slots_mutex));
    /* Legacy results expire at the next encode, never prepared results. */
    if (writer->legacy) {
        writer->legacy->in_use = 0;
        writer->legacy = NULL;
    }
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i) {
        if (!writer->slots[i].in_use) {
            batch = &writer->slots[i];
            batch->in_use = 1;
            break;
        }
    }
    mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
    if (!batch) errno = EAGAIN;
    return batch;
}

int swbwa_gzip_output_release(swbwa_gzip_output *writer, swbwa_gzip_output_batch *batch)
{
    unsigned i;
    int result = -1;
    if (!batch) return 0;
    if (!writer) { errno = EINVAL; return -1; }
    mutex_check(pthread_mutex_lock(&writer->slots_mutex));
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i) {
        if (batch == &writer->slots[i] && batch->in_use) {
            batch->in_use = 0;
            if (writer->legacy == batch) writer->legacy = NULL;
            result = 0;
            break;
        }
    }
    mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
    if (result) errno = EINVAL;
    return result;
}

int swbwa_gzip_output_pending(swbwa_gzip_output *writer)
{
    unsigned i;
    int count = 0;
    mutex_check(pthread_mutex_lock(&writer->slots_mutex));
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i)
        count += writer->slots[i].in_use && &writer->slots[i] != writer->legacy;
    mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
    return count;
}

int swbwa_gzip_output_batch_data(swbwa_gzip_output *writer,
                                 swbwa_gzip_output_batch *batch, int reads,
                                 const struct iovec **vectors, unsigned *count,
                                 uint64_t *raw_bytes, uint64_t *gzip_bytes)
{
    unsigned i;
    int result = -1;
    if (!writer || !vectors || !count || !raw_bytes || !gzip_bytes) {
        errno = EINVAL;
        return -1;
    }
    mutex_check(pthread_mutex_lock(&writer->slots_mutex));
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i) {
        if (batch == &writer->slots[i] && batch->in_use && batch->reads == reads) {
            *vectors = batch->vectors;
            *count = batch->count;
            *raw_bytes = batch->raw_bytes;
            *gzip_bytes = batch->gzip_bytes;
            result = 0;
            break;
        }
    }
    mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
    if (result) errno = EINVAL;
    return result;
}

enum { GZIP_WORKSPACE_BYTES = 160 * 1024 };

static double monotonic_seconds(void)
{
    struct timespec time;
    clock_gettime(CLOCK_MONOTONIC, &time);
    return time.tv_sec + time.tv_nsec * 1.e-9;
}

static int environment_integer(const char *name, long minimum, long maximum, long *value)
{
    const char *text = getenv(name);
    char *end;
    long parsed;
    if (!text || !*text) return 0;
    errno = 0;
    parsed = strtol(text, &end, 10);
    if (errno || *end || parsed < minimum || parsed > maximum) {
        fprintf(stderr, "[E::gzip output] %s must be in %ld..%ld\n", name, minimum, maximum);
        errno = EINVAL;
        return -1;
    }
    *value = parsed;
    return 0;
}

swbwa_gzip_output *swbwa_gzip_output_create(unsigned workers)
{
    long level = 1, kib = 0;
    swbwa_gzip_output *writer;
    if (!workers || workers > SWBWA_CPE_COUNT) { errno = EINVAL; return NULL; }
    if (environment_integer("SWBWA_GZIP_LEVEL", 1, 9, &level) ||
        environment_integer("SWBWA_GZIP_BLOCK_KIB", 0, 4096, &kib)) return NULL;
    writer = calloc(1, sizeof(*writer));
    if (!writer) return NULL;
    writer->workers = workers;
    writer->level = (int)level;
    writer->block_bytes = (size_t)kib * 1024;
    writer->workspace_storage = malloc((size_t)workers * GZIP_WORKSPACE_BYTES + 63);
    if (!writer->workspace_storage) { free(writer); return NULL; }
    writer->workspaces = (unsigned char *)(((uintptr_t)writer->workspace_storage + 63) & ~(uintptr_t)63);
    mutex_check(pthread_mutex_init(&writer->encoder_mutex, NULL));
    mutex_check(pthread_mutex_init(&writer->slots_mutex, NULL));
    return writer;
}

static void measure_records(void *opaque, int worker, int workers)
{
    swbwa_gzip_output *writer = opaque;
    size_t begin = (size_t)writer->reads * worker / workers;
    size_t end = (size_t)writer->reads * (worker + 1) / workers;
    size_t total = 0, i;
    for (i = begin; i < end; ++i) {
        size_t bytes;
        if (!writer->seqs[i].sam || (writer->lengths && writer->lengths[i] < 0)) {
            writer->errors[worker] = EINVAL;
            return;
        }
        bytes = writer->lengths ? (size_t)writer->lengths[i] : strlen(writer->seqs[i].sam);
        if (bytes > SIZE_MAX - total) { writer->errors[worker] = EOVERFLOW; return; }
        writer->spans[i].data = (const unsigned char *)writer->seqs[i].sam;
        writer->spans[i].size = bytes;
        total += bytes;
    }
    writer->sums[worker] = total;
}

static void prefix_records(void *opaque, int worker, int workers)
{
    swbwa_gzip_output *writer = opaque;
    size_t begin = (size_t)writer->reads * worker / workers;
    size_t end = (size_t)writer->reads * (worker + 1) / workers;
    size_t position = writer->offsets[worker], i;
    for (i = begin; i < end; ++i) {
        writer->prefix[i] = position;
        position += writer->spans[i].size;
    }
}

static int reserve_records(swbwa_gzip_output *writer, size_t records)
{
    swbwa_codec_span *spans;
    size_t *prefix;
    if (records <= writer->record_capacity) return 0;
    if (records >= SIZE_MAX / sizeof(*prefix) || records > SIZE_MAX / sizeof(*spans)) {
        errno = EOVERFLOW;
        return -1;
    }
    spans = realloc(writer->spans, records * sizeof(*spans));
    if (!spans) return -1;
    writer->spans = spans;
    prefix = realloc(writer->prefix, (records + 1) * sizeof(*prefix));
    if (!prefix) return -1;
    writer->prefix = prefix;
    writer->record_capacity = records;
    return 0;
}

static size_t record_at(const swbwa_gzip_output *writer, size_t bytes)
{
    size_t left = 0, right = (size_t)writer->reads;
    while (left < right) {
        size_t middle = left + (right - left) / 2;
        if (writer->prefix[middle + 1] <= bytes) left = middle + 1;
        else right = middle;
    }
    return left;
}

static int encode(swbwa_gzip_output *writer, swbwa_gzip_output_batch *output,
                   size_t bytes, double start)
{
    size_t block = writer->block_bytes, stride, needed, position;
    unsigned jobs, i;
    swbwa_codec_batch batch;
    if (!block) {
        block = bytes / ((size_t)writer->workers * 2);
        if (block < 16384) block = 16384;
        if (block > 262144) block = 262144;
        block = (block + 4095) & ~(size_t)4095;
    }
    if (bytes / block >= UINT_MAX) { errno = EOVERFLOW; return -1; }
    jobs = bytes ? (unsigned)(1 + (bytes - 1) / block) : 1;
    /* Conservative bound for arbitrary zlib memory settings, plus alignment. */
    stride = (block + (block >> 3) + (block >> 6) + 256 + 63) & ~(size_t)63;
    if (jobs > SIZE_MAX / stride) { errno = EOVERFLOW; return -1; }
    needed = (size_t)jobs * stride;
    if (jobs > writer->job_capacity) {
        swbwa_codec_job *new_jobs = realloc(writer->jobs, (size_t)jobs * sizeof(*new_jobs));
        if (!new_jobs) return -1;
        writer->jobs = new_jobs;
        writer->job_capacity = jobs;
    }
    if (jobs > output->vector_capacity) {
        struct iovec *new_vectors = realloc(output->vectors, (size_t)jobs * sizeof(*new_vectors));
        if (!new_vectors) return -1;
        output->vectors = new_vectors;
        output->vector_capacity = jobs;
    }
    if (needed > output->slab_capacity) {
        unsigned char *slab = realloc(output->slab, needed);
        if (!slab) return -1;
        output->slab = slab;
        output->slab_capacity = needed;
    }
    for (i = 0, position = 0; i < jobs; ++i) {
        swbwa_codec_job *job = &writer->jobs[i];
        size_t length = bytes - position < block ? bytes - position : block;
        size_t first = bytes ? record_at(writer, position) : 0;
        size_t last = length ? record_at(writer, position + length - 1) : first;
        memset(job, 0, sizeof(*job));
        job->operation = SWBWA_CODEC_GZIP;
        job->level = writer->level;
        job->input_size = length;
        job->spans = writer->spans + first;
        job->span_count = length ? (unsigned)(last - first + 1) : 0;
        job->span_skip = bytes ? position - writer->prefix[first] : 0;
        job->output = output->slab + (size_t)i * stride;
        job->output_capacity = stride;
        position += length;
    }
    writer->prepare_seconds += monotonic_seconds() - start;
    memset(&batch, 0, sizeof(batch));
    batch.jobs = writer->jobs;
    batch.count = jobs;
    batch.workers = writer->workers;
    batch.workspaces = writer->workspaces;
    batch.workspace_bytes = GZIP_WORKSPACE_BYTES;
    if ((jobs < writer->workers ? jobs : writer->workers) > writer->max_active_workers)
        writer->max_active_workers = jobs < writer->workers ? jobs : writer->workers;
    swbwa_codec_launch(&batch);
    writer->cpe_seconds += batch.run_s;
    writer->wait_seconds += batch.wait_s;
    output->gzip_bytes = 0;
    for (i = 0; i < jobs; ++i) {
        swbwa_codec_job *job = &writer->jobs[i];
        if (job->status != 0 || !job->output_size) {
            fprintf(stderr, "[E::gzip output] member=%u status=%d input=%zu output=%zu\n",
                    i, job->status, job->input_size, job->output_size);
            errno = EIO;
            return -1;
        }
        output->vectors[i].iov_base = job->output;
        output->vectors[i].iov_len = job->output_size;
        output->gzip_bytes += job->output_size;
        writer->ldm_members += job->state_in_ldm != 0;
    }
    output->count = jobs;
    output->reads = writer->reads;
    output->raw_bytes = bytes;
    ++writer->batches;
    writer->members += jobs;
    writer->input_bytes += bytes;
    writer->output_bytes += output->gzip_bytes;
    return 0;
}

static int encode_records(swbwa_gzip_output *writer, swbwa_gzip_output_batch *output,
                           const bseq1_t *seqs, const int *lengths, int reads)
{
    double start = monotonic_seconds();
    size_t bytes = 0;
    int i;
    if (reserve_records(writer, reads ? (size_t)reads : 1)) return -1;
    writer->seqs = seqs;
    writer->lengths = lengths;
    writer->reads = reads;
    memset(writer->errors, 0, sizeof(writer->errors));
    memset(writer->sums, 0, sizeof(writer->sums));
    swbwa_host_workers_run_ready(measure_records, writer);
    for (i = 0; i < SWBWA_HOST_MPE_THREADS; ++i) {
        if (writer->errors[i] || writer->sums[i] > SIZE_MAX - bytes) {
            errno = writer->errors[i] ? writer->errors[i] : EOVERFLOW;
            return -1;
        }
        writer->offsets[i] = bytes;
        bytes += writer->sums[i];
    }
    swbwa_host_workers_run_ready(prefix_records, writer);
    writer->prefix[reads] = bytes;
    return encode(writer, output, bytes, start);
}

static int prepare_records(swbwa_gzip_output *writer, const bseq1_t *seqs,
                             const int *lengths, int reads, swbwa_gzip_output_batch **batch)
{
    swbwa_gzip_output_batch *output;
    int result = -1, saved;
    output = acquire_batch(writer);
    if (output) {
        result = encode_records(writer, output, seqs, lengths, reads);
        if (!result) *batch = output;
        else {
            saved = errno;
            swbwa_gzip_output_release(writer, output);
            errno = saved;
        }
    }
    return result;
}

int swbwa_gzip_output_prepare(swbwa_gzip_output *writer, const bseq1_t *seqs,
                              const int *lengths, int reads, swbwa_gzip_output_batch **batch)
{
    int result, saved;
    if (batch) *batch = NULL;
    if (!writer || !batch || reads < 0 || (reads && !seqs)) { errno = EINVAL; return -1; }
    if (encoder_lock(writer)) return -1;
    result = prepare_records(writer, seqs, lengths, reads, batch);
    saved = errno;
    mutex_check(pthread_mutex_unlock(&writer->encoder_mutex));
    errno = saved;
    return result;
}

int swbwa_gzip_output_encode(swbwa_gzip_output *writer, const bseq1_t *seqs,
                             const int *lengths, int reads,
                             const struct iovec **vectors, unsigned *count,
                             uint64_t *raw_bytes, uint64_t *gzip_bytes)
{
    swbwa_gzip_output_batch *output;
    int result, saved;
    if (!writer || reads < 0 || (reads && !seqs) || !vectors || !count || !raw_bytes || !gzip_bytes) {
        errno = EINVAL;
        return -1;
    }
    if (encoder_lock(writer)) return -1;
    result = prepare_records(writer, seqs, lengths, reads, &output);
    if (!result) {
        mutex_check(pthread_mutex_lock(&writer->slots_mutex));
        writer->legacy = output;
        mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
        result = swbwa_gzip_output_batch_data(writer, output, reads, vectors, count, raw_bytes, gzip_bytes);
    }
    saved = errno;
    mutex_check(pthread_mutex_unlock(&writer->encoder_mutex));
    errno = saved;
    return result;
}

int swbwa_gzip_output_bytes(swbwa_gzip_output *writer, const void *data, size_t bytes,
                            const struct iovec **vectors, unsigned *count,
                            uint64_t *gzip_bytes)
{
    double start = monotonic_seconds();
    swbwa_gzip_output_batch *output;
    int result = -1, saved;
    if (!writer || (bytes && !data) || !vectors || !count || !gzip_bytes) {
        errno = EINVAL;
        return -1;
    }
    if (encoder_lock(writer)) return -1;
    output = acquire_batch(writer);
    if (output) {
        if (!reserve_records(writer, 1)) {
            writer->reads = 1;
            writer->spans[0].data = data;
            writer->spans[0].size = bytes;
            writer->prefix[0] = 0;
            writer->prefix[1] = bytes;
            result = encode(writer, output, bytes, start);
        }
        saved = errno;
        if (result) swbwa_gzip_output_release(writer, output);
        else {
            mutex_check(pthread_mutex_lock(&writer->slots_mutex));
            writer->legacy = output;
            mutex_check(pthread_mutex_unlock(&writer->slots_mutex));
            *vectors = output->vectors;
            *count = output->count;
            *gzip_bytes = output->gzip_bytes;
        }
        errno = saved;
    }
    saved = errno;
    mutex_check(pthread_mutex_unlock(&writer->encoder_mutex));
    errno = saved;
    return result;
}

void swbwa_gzip_output_report(const swbwa_gzip_output *writer, FILE *file, int rank, int ranks)
{
    fprintf(file, "[SWBWA gzip output rank %06d/%06d] workers=%u level=%d block_kib=%zu"
            " batches=%" PRIu64 " members=%" PRIu64 " raw_bytes=%" PRIu64
            " gzip_bytes=%" PRIu64 " ratio=%.3f prepare_s=%.6f cpe_s=%.6f wait_s=%.6f"
            " ldm_members=%" PRIu64 " max_active_workers=%u mpe_sam_copy_bytes=0\n", rank, ranks,
            writer->workers, writer->level, writer->block_bytes / 1024,
            writer->batches, writer->members, writer->input_bytes, writer->output_bytes,
            writer->output_bytes ? (double)writer->input_bytes / writer->output_bytes : 0.,
            writer->prepare_seconds, writer->cpe_seconds, writer->wait_seconds, writer->ldm_members,
            writer->max_active_workers);
}

void swbwa_gzip_output_destroy(swbwa_gzip_output *writer)
{
    unsigned i;
    if (!writer) return;
    free(writer->spans);
    free(writer->prefix);
    free(writer->jobs);
    for (i = 0; i < SWBWA_PIPELINE_BUFFER_COUNT; ++i) {
        free(writer->slots[i].vectors);
        free(writer->slots[i].slab);
    }
    mutex_check(pthread_mutex_destroy(&writer->encoder_mutex));
    mutex_check(pthread_mutex_destroy(&writer->slots_mutex));
    free(writer->workspace_storage);
    free(writer);
}
