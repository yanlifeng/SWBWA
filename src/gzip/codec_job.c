#include "codec_job.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#if defined(SWBWA_CODEC_CPE)
#include <slave.h>
#endif
#ifdef SWBWA_CODEC_RAPID
extern void swbwa_rapid_discover(swbwa_codec_job *job);
#endif

typedef struct {
    unsigned char *local;
    size_t used, capacity;
    unsigned char *main;
    size_t main_used, main_capacity;
} gzip_workspace;

static voidpf gzip_allocate(voidpf opaque, uInt count, uInt size)
{
    gzip_workspace *workspace = opaque;
    size_t bytes, aligned;
    void *pointer;
    if (size && count > SIZE_MAX / size) return NULL;
    bytes = (size_t)count * size;
    if (bytes > SIZE_MAX - 63) return NULL;
    aligned = (bytes + 63) & ~(size_t)63;
    if (aligned <= workspace->capacity - workspace->used) {
        pointer = workspace->local + workspace->used;
        workspace->used += aligned;
        memset(pointer, 0, bytes);
        return pointer;
    }
    if (aligned <= workspace->main_capacity - workspace->main_used) {
        pointer = workspace->main + workspace->main_used;
        workspace->main_used += aligned;
        memset(pointer, 0, bytes);
        return pointer;
    }
    return calloc(count, size);
}

static void gzip_release(voidpf opaque, voidpf pointer)
{
    gzip_workspace *workspace = opaque;
    uintptr_t address = (uintptr_t)pointer, local = (uintptr_t)workspace->local;
    uintptr_t arena = (uintptr_t)workspace->main;
    if (workspace->main && address >= arena && address - arena < workspace->main_capacity)
        return;
    if (!workspace->local || address < local || address - local >= workspace->capacity)
        free(pointer);
}

/* Reuse one zlib state for every job assigned to this CPE in the dispatch.
 * Small state/hash/pending arrays fit LDM; the 64 KiB window/prev arrays fall
 * back to main memory. No allocator state survives a CPE launch. */
void swbwa_codec_gzip_worker(swbwa_codec_batch *batch, unsigned worker)
{
    z_stream stream;
    gzip_workspace workspace = {0};
    unsigned index;
    int initialized = 0, level = -1;
#if defined(SWBWA_CODEC_CPE)
    void *backing = NULL;
#endif
    if (worker >= batch->workers || worker >= batch->count) return;
    if (batch->workspaces) {
        workspace.main = batch->workspaces + (size_t)worker * batch->workspace_bytes;
        workspace.main_capacity = batch->workspace_bytes;
    }
    memset(&stream, 0, sizeof(stream));
#if defined(SWBWA_CODEC_CPE)
    backing = ldm_malloc(32768 + 63);
    if (backing) {
        workspace.local = (unsigned char *)(((uintptr_t)backing + 63) & ~(uintptr_t)63);
        workspace.capacity = 32768;
    }
#endif
    stream.zalloc = gzip_allocate;
    stream.zfree = gzip_release;
    stream.opaque = &workspace;
    for (index = worker; index < batch->count; index += batch->workers) {
        swbwa_codec_job *job = &batch->jobs[index];
        size_t remaining = job->input_size, skip = job->span_skip;
        unsigned span;
        int result;
        job->status = Z_STREAM_ERROR;
        job->output_size = 0;
        if (job->operation != SWBWA_CODEC_GZIP || !job->output ||
            job->output_capacity > UINT_MAX || job->input_size > UINT_MAX ||
            (remaining && (!job->spans || !job->span_count))) continue;
        if (initialized && level != job->level) {
            deflateEnd(&stream);
            initialized = 0;
            workspace.used = 0;
            workspace.main_used = 0;
        }
        if (!initialized) {
            result = deflateInit2(&stream, job->level, Z_DEFLATED, 15 + 16,
                                  5, Z_DEFAULT_STRATEGY);
            if (result != Z_OK) {
                job->status = result;
                workspace.used = workspace.main_used = 0;
                continue;
            }
            initialized = 1;
            level = job->level;
        } else if ((result = deflateReset(&stream)) != Z_OK) {
            job->status = result;
            continue;
        }
        stream.next_out = job->output;
        stream.avail_out = (uInt)job->output_capacity;
        result = Z_OK;
        for (span = 0; span < job->span_count && remaining; ++span) {
            const swbwa_codec_span *part = &job->spans[span];
            size_t bytes;
            if (skip > part->size || (part->size && !part->data)) {
                result = Z_STREAM_ERROR;
                break;
            }
            bytes = part->size - skip;
            if (bytes > remaining) bytes = remaining;
            if (!bytes) { skip = 0; continue; }
            stream.next_in = (Bytef *)(part->data + skip);
            stream.avail_in = (uInt)bytes;
            do {
                result = deflate(&stream, Z_NO_FLUSH);
            } while (result == Z_OK && stream.avail_in && stream.avail_out);
            if (result != Z_OK || stream.avail_in || !stream.avail_out) {
                if (result == Z_OK) result = Z_BUF_ERROR;
                break;
            }
            remaining -= bytes;
            skip = 0;
        }
        if (result == Z_OK && remaining) result = Z_STREAM_ERROR;
        if (result == Z_OK) result = deflate(&stream, Z_FINISH);
        job->output_size = job->output_capacity - stream.avail_out;
        job->state_in_ldm = workspace.used != 0;
        job->state_bytes = workspace.used;
        job->status = result == Z_STREAM_END ? Z_OK :
                      result == Z_OK ? Z_BUF_ERROR : result;
    }
    if (initialized) deflateEnd(&stream);
#if defined(SWBWA_CODEC_CPE)
    if (backing) ldm_free(backing, 32768 + 63);
#endif
}

void swbwa_codec_execute(swbwa_codec_job *job)
{
    z_stream stream;
    int result;
    memset(&stream, 0, sizeof(stream));
    job->status = Z_STREAM_ERROR;
    if (job->operation == SWBWA_CODEC_RESOLVE) {
        size_t i;
        if (!job->symbols || !job->output || job->dictionary_size > 32768 ||
            job->output_size > job->output_capacity || job->output_size > UINT_MAX) return;
        for (i = 0; i < job->output_size;) {
            size_t page = i / 4096;
            size_t end = job->output_size - i < 4096 ? job->output_size : i + 4096;
            uint32_t lines = 0;
            for (; i < end; ++i) {
                unsigned symbol = job->symbols[i];
                if (symbol >= 256) {
                    if (symbol < 65536U - job->dictionary_size || !job->dictionary) {
                        job->status = Z_DATA_ERROR;
                        return;
                    }
                    symbol = job->dictionary[symbol - (65536U - job->dictionary_size)];
                }
                job->output[i] = (unsigned char)symbol;
                lines += symbol == '\n';
            }
            if (job->line_counts) job->line_counts[page] = lines;
        }
        job->crc = (uint32_t)crc32(0, job->output, (uInt)job->output_size);
        job->status = Z_OK;
        return;
    }
#ifdef SWBWA_CODEC_RAPID
    if (job->operation == SWBWA_CODEC_DISCOVER) {
        swbwa_rapid_discover(job);
        return;
    }
#endif
    job->output_size = 0;
    if (job->operation != SWBWA_CODEC_DEFLATE || !job->input || !job->output || job->input_size > UINT_MAX ||
        job->output_capacity > UINT_MAX || job->dictionary_size > 32768 ||
        (job->dictionary_size && !job->dictionary))
        return;
    stream.next_in = (Bytef *)job->input;
    stream.avail_in = (uInt)job->input_size;
    stream.next_out = job->output;
    stream.avail_out = (uInt)job->output_capacity;
    result = deflateInit2(&stream, job->level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
    if (result != Z_OK) { job->status = result; return; }
    if (job->dictionary_size)
        result = deflateSetDictionary(&stream, job->dictionary, job->dictionary_size);
    if (result == Z_OK) result = deflate(&stream, Z_SYNC_FLUSH);
    job->output_size = job->output_capacity - stream.avail_out;
    job->status = result == Z_OK && !stream.avail_in && stream.avail_out ? Z_OK : Z_BUF_ERROR;
    if (job->status == Z_OK)
        job->crc = (uint32_t)crc32(0, job->input, (uInt)job->input_size);
    /* An intentionally unfinished raw stream ends with Z_SYNC_FLUSH.
     * The host appends one final block after concatenating all jobs. */
    deflateEnd(&stream);
}
