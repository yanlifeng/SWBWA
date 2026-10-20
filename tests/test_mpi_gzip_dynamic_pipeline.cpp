#include "swbwa_mpi_gzip.h"
#include "swbwa_pipeline.h"
#include "swbwa_host_workers.h"
#include "codec_job.h"
#include "mpi_decoder.h"
extern "C" {
#include "swbwa_mpi.h"
#include "swbwa_output.h"
void kt_pipeline_single(int, void *(*)(void *, int, void *), void *, int);
void kt_pipeline_queue(int, void *(*)(void *, int, void *), void *, int);
}
#include <cassert>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <pthread.h>
#include <string>
#include <thread>
#include <vector>

static std::thread::id owner;
struct ThreadProbe {
    std::thread::id expected;
    unsigned calls = 0;
    explicit ThreadProbe(std::thread::id thread) : expected(thread) {}
};
static void probe_thread(void *opaque)
{
    ThreadProbe &p = *static_cast<ThreadProbe *>(opaque);
    assert(std::this_thread::get_id() == p.expected);
    ++p.calls;
}

static std::atomic<bool> capture_queue_key{false}, queue_key_captured{false};
static pthread_key_t queue_key;
static std::atomic<unsigned> key_sets{0}, key_clears{0}, clear_direct_calls{0};

extern "C" int __real_pthread_key_create(pthread_key_t *, void (*)(void *));
extern "C" int __wrap_pthread_key_create(pthread_key_t *key, void (*destructor)(void *))
{
    const int error = __real_pthread_key_create(key, destructor);
    if (capture_queue_key.load() && std::this_thread::get_id() == owner) {
        assert(!error && !destructor && !queue_key_captured.load());
        queue_key = *key;
        queue_key_captured.store(true);
    }
    return error;
}

extern "C" int __real_pthread_setspecific(pthread_key_t, const void *);
extern "C" int __wrap_pthread_setspecific(pthread_key_t key, const void *value)
{
    const bool watched = queue_key_captured.load() && key == queue_key;
    if (watched) {
        assert(std::this_thread::get_id() != owner);
        assert(value ? pthread_getspecific(key) == NULL : pthread_getspecific(key) != NULL);
    }
    const int error = __real_pthread_setspecific(key, value);
    if (watched) {
        assert(!error && pthread_getspecific(key) == value);
        if (value) ++key_sets;
        else {
            ++key_clears;
            // Check explicit cleanup before thread exit can discard its value.
            ThreadProbe p(std::this_thread::get_id());
            swbwa_pipeline_call_on_processor(probe_thread, &p);
            assert(p.calls == 1);
            ++clear_direct_calls;
        }
    }
    return error;
}

struct HandoffFixture {
    ThreadProbe processor{owner};
    unsigned total, produced = 0, processed = 0, written = 0, flushed = 0;
    unsigned tokens[3] = {0, 1, 2};
    explicit HandoffFixture(unsigned count) : total(count) {}
};

static void *handoff_process(void *opaque, int step, void *item)
{
    HandoffFixture &f = *static_cast<HandoffFixture *>(opaque);
    if (step == 0) {
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        return f.produced == f.total ? NULL : &f.tokens[f.produced++];
    }
    if (step == 1) {
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        assert(*static_cast<unsigned *>(item) == f.processed++);
        return item;
    }
    if (step == 2) {
        ThreadProbe writer(std::this_thread::get_id());
        swbwa_pipeline_call_on_processor(probe_thread, &writer);
        assert(writer.calls == 1 && *static_cast<unsigned *>(item) == f.written++);
        return NULL;
    }
    assert(step == 3);
    ++f.flushed;
    return NULL;
}

static void check_handoff_reuse()
{
    for (unsigned count : {0U, 3U}) {
        HandoffFixture f(count);
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        kt_pipeline_queue(3, handoff_process, &f, 3);
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        assert(f.produced == count && f.processed == count && f.written == count && f.flushed == 1);
        assert(f.processor.calls == 2 + 2 * (count + 1) + count);
    }
}

static std::string stress;
static int ranks;
static unsigned ready_limit;
static struct {
    std::atomic<unsigned> pending{0}, overlap_refills{0}, codec_with_pending{0};
    uint64_t decode_calls = 0, decodes_after_overlap = 0;
    uint64_t window_decode_calls = 0, overlap_streak = 0, max_overlap_streak = 0;
    uint64_t read_chunks = 0, processed_chunks = 0, written_chunks = 0;
    int64_t last_read = -1, last_processed = -1, last_written = -1, last_window = -1;
} observed;

static void check_decode_drain()
{
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
    if (ranks > 1) assert(observed.pending.load() == 0);
#endif
}

extern "C" int __real_swbwa_mpi_decoder_next(swbwa_mpi_decoder *, swbwa_mpi_decoded_batch *);
extern "C" int __wrap_swbwa_mpi_decoder_next(swbwa_mpi_decoder *d, swbwa_mpi_decoded_batch *b)
{
    check_decode_drain();
    ++observed.decode_calls;
    if (observed.max_overlap_streak >= 2) ++observed.decodes_after_overlap;
    return __real_swbwa_mpi_decoder_next(d, b);
}

extern "C" void __real_swbwa_codec_launch(swbwa_codec_batch *);
static void launch(void *p)
{
    assert(std::this_thread::get_id() == owner);
    check_decode_drain();
    if (observed.pending.load()) ++observed.codec_with_pending;
    __real_swbwa_codec_launch(static_cast<swbwa_codec_batch *>(p));
}
extern "C" void __wrap_swbwa_codec_launch(swbwa_codec_batch *b)
{
    swbwa_pipeline_call_on_processor(launch, b);
}
extern "C" void __real_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *);
extern "C" void __wrap_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *c)
{
    __real_swbwa_mpi_decoder_config_init(c);
    c->chunk_bytes = 4096; c->jobs_per_rank = 1; c->initial_output_bytes = 65536;
    c->decoded_bytes = 8U << 20; c->memory_bytes = 64U << 20; c->lookahead_bytes = 1U << 20;
}
struct Batch {
    swbwa_mpi_gzip_batch input = {};
    std::vector<std::string> text;
    std::vector<bseq1_t> seqs;
    ~Batch()
    {
        std::free(input.data[0]); std::free(input.data[1]); std::free(input.partitions);
    }
};
struct Fixture { swbwa_mpi_gzip_reader *reader; int mates, rank; };

static void wait_for_overlap()
{
    // The first writer releases only after the reader proves the intended overlap.
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    while (stress == "decoded-refills" ? observed.overlap_refills.load() < 2 :
                                       observed.codec_with_pending.load() == 0) {
        if (std::chrono::steady_clock::now() >= deadline) {
            std::fprintf(stderr, "pipeline overlap gate timed out: %s\n", stress.c_str());
            MPI_Abort(MPI_COMM_WORLD, 1);
            std::abort();
        }
        assert(!swbwa_mpi_progress());
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

static void *process(void *context, int step, void *item)
{
    Fixture &f = *static_cast<Fixture *>(context);
    Batch *b = static_cast<Batch *>(item);
    if (step == 0) {
        b = new Batch;
        const int status = swbwa_mpi_gzip_next(f.reader, &b->input);
        assert(status >= 0);
        if (!status) { delete b; return NULL; }
        assert(b->input.chunk_id > observed.last_read);
        observed.last_read = b->input.chunk_id;
        if (b->input.records) {
#if SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
            const int64_t window = b->input.chunk_id / (int64_t(ranks) * ready_limit);
            if (ranks > 1 && window != observed.last_window) {
                if (stress == "window-growth")
                    assert(!swbwa_mpi_gzip_report(f.reader, stderr));
                if (observed.last_window >= 0 &&
                    observed.decode_calls == observed.window_decode_calls && observed.pending.load()) {
                    ++observed.overlap_refills;
                    ++observed.overlap_streak;
                    if (observed.overlap_streak > observed.max_overlap_streak)
                        observed.max_overlap_streak = observed.overlap_streak;
                } else observed.overlap_streak = 0;
                observed.last_window = window;
                observed.window_decode_calls = observed.decode_calls;
            }
#endif
            ++observed.pending;
            ++observed.read_chunks;
        }
        return b;
    }
    if (step == 1) {
        assert(std::this_thread::get_id() == owner);
        assert(b->input.chunk_id > observed.last_processed);
        observed.last_processed = b->input.chunk_id;
        if (b->input.records) ++observed.processed_chunks;
        if (!stress.empty() && f.rank == 0 &&
            b->input.chunk_id % (int64_t(ranks) * ready_limit) == 0)
            std::this_thread::sleep_for(std::chrono::milliseconds(8));
        if (b->input.records && (b->input.chunk_id % 7 == 0 || f.rank == 1))
            std::this_thread::sleep_for(std::chrono::milliseconds(3));
        size_t offsets[2] = {0, 0};
        int64_t ordinal = b->input.first_record;
        for (unsigned c = 0; c < b->input.partition_count; ++c) {
            const swbwa_fastq_partition_t &p = b->input.partitions[c];
            assert(ordinal == b->input.first_record + p.first_record);
            for (int m = 0; m < f.mates; ++m) assert(size_t(p.begin[m]) == offsets[m]);
            for (int64_t i = 0; i < p.records; ++i, ++ordinal) {
                std::string name;
                for (int m = 0; m < f.mates; ++m) {
                    swbwa_fastq_record_view v;
                    assert(swbwa_fastq_record_next(b->input.data[m], size_t(p.end[m]), &offsets[m], &v) == 1);
                    const std::string current(v.name, v.name_size);
                    if (m) assert(name == current); else name = current;
                    b->text.push_back(current + "\t4\t*\t0\t0\t*\t*\t0\t0\t" +
                        std::string(v.seq, v.length) + "\t" + std::string(v.qual, v.length) +
                        "\tGI:i:" + std::to_string(ordinal * f.mates + m) + "\n");
                }
            }
            for (int m = 0; m < f.mates; ++m) assert(offsets[m] == size_t(p.end[m]));
        }
        assert(ordinal == b->input.first_record + b->input.records);
        b->seqs.resize(b->text.size());
        for (size_t i = 0; i < b->text.size(); ++i) b->seqs[i].sam = const_cast<char *>(b->text[i].c_str());
        return b;
    }
    if (step == 2) {
        if (b->input.records) {
            assert(b->input.chunk_id > observed.last_written);
            observed.last_written = b->input.chunk_id;
            ++observed.written_chunks;
            if (b->input.chunk_id == 0 && (stress == "decoded-refills" || stress == "single-direct"))
                wait_for_overlap();
            if (!stress.empty() && f.rank == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(4));
            if (b->input.chunk_id % 5 == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            assert(swbwa_output_write_chunk(b->input.chunk_id, b->input.first_record,
                b->input.first_record + b->input.records, b->seqs.data(), NULL,
                int(b->seqs.size()), f.mates == 2) == 0);
            // Decrement before complete(), which can immediately wake a decoding reader.
            assert(observed.pending.fetch_sub(1) > 0);
            assert(!swbwa_mpi_gzip_complete(f.reader, b->input.chunk_id));
        }
        delete b;
        return NULL;
    }
    assert(step == 3);
    assert(!swbwa_output_flush());
    return NULL;
}

static void reject_pending_close(swbwa_mpi_gzip_reader *reader, int rank)
{
    std::vector<std::unique_ptr<Batch>> held;
    for (;;) {
        std::unique_ptr<Batch> b(new Batch);
        const int status = swbwa_mpi_gzip_next(reader, &b->input);
        assert(status >= 0);
        if (!status) break;
        assert(b->input.records > 0);
        ++observed.pending;
        held.push_back(std::move(b));
    }
    assert(swbwa_mpi_gzip_chunk_count(reader) == 1);
    assert(held.size() == size_t(rank == 0));
    std::fprintf(stderr, "[test pending close rank=%d] validated_EOF=1 pending=%zu\n",
                 rank, held.size());
    std::fflush(stderr);
    swbwa_mpi_gzip_close(reader);
    std::fprintf(stderr, "pending close unexpectedly returned\n");
    std::abort();
}

int main(int argc, char **argv)
{
    assert(argc == 5);
    assert(!swbwa_mpi_init(&argc, &argv));
    ranks = swbwa_mpi_size();
    const char *profile = std::getenv("PIPELINE_TEST_STRESS");
    stress = profile ? profile : "";
    assert(stress.empty() || stress == "decoded-refills" || stress == "window-reuse" ||
           stress == "single-direct" || stress == "window-growth" || stress == "pending-close");
    const char *ready = std::getenv("SWBWA_GZIP_READY_CHUNKS");
    ready_limit = ready ? unsigned(std::strtoul(ready, NULL, 10)) : 4;
    assert(ready_limit >= 1 && ready_limit <= 16);
    owner = std::this_thread::get_id();
    ThreadProbe initial(owner);
    capture_queue_key.store(true);
    swbwa_pipeline_call_on_processor(probe_thread, &initial);
    capture_queue_key.store(false);
    assert(queue_key_captured.load() && initial.calls == 1);
    assert(!swbwa_host_workers_init());
    swbwa_codec_init();
    Fixture f = {NULL, std::strcmp(argv[2], "-") ? 2 : 1, swbwa_mpi_rank()};
    f.reader = swbwa_mpi_gzip_open(argv[1], f.mates == 2 ? argv[2] : NULL,
                                 stress == "window-growth" ? 1 : stress.empty() ? 2048 : 1024, 64, 1);
    assert(f.reader);
    if (stress == "pending-close") reject_pending_close(f.reader, f.rank);
    assert(!swbwa_output_open_online(argv[3], 1));
    assert(!swbwa_mpi_progress_thread_start(0));
    if (std::strcmp(argv[4], "queue") == 0) kt_pipeline_queue(3, process, &f, 3);
    else kt_pipeline_single(1, process, &f, 3);
    assert(observed.pending.load() == 0);
    assert(observed.read_chunks == observed.processed_chunks &&
           observed.read_chunks == observed.written_chunks);
    if (stress == "decoded-refills" && !f.rank) {
        assert(observed.max_overlap_streak >= 2);
        assert(observed.decodes_after_overlap > 0);
    }
    if (stress == "single-direct") assert(observed.codec_with_pending.load() > 0);
    std::fprintf(stderr,
        "[test pipeline rank=%d] chunks=%llu decode_calls=%llu overlap_refills=%u "
        "max_overlap_streak=%llu decodes_after_overlap=%llu codec_with_pending=%u\n",
        f.rank, (unsigned long long)observed.read_chunks, (unsigned long long)observed.decode_calls,
        observed.overlap_refills.load(), (unsigned long long)observed.max_overlap_streak,
        (unsigned long long)observed.decodes_after_overlap, observed.codec_with_pending.load());
    assert(!swbwa_output_ordered_finish(swbwa_mpi_gzip_chunk_count(f.reader)));
    assert(!swbwa_mpi_gzip_report(f.reader, stderr));
    swbwa_mpi_gzip_close(f.reader);
    assert(!swbwa_output_close());
    assert(!swbwa_mpi_progress_thread_stop());
    check_handoff_reuse();
    const unsigned pipelines = 2 + unsigned(std::strcmp(argv[4], "queue") == 0);
    assert(key_sets.load() == pipelines && key_clears.load() == pipelines &&
           clear_direct_calls.load() == pipelines);
    std::fprintf(stderr,
        "[test pipeline keys rank=%d] pipelines=%u sets=%u clears=%u clear_direct_calls=%u\n",
        f.rank, pipelines, key_sets.load(), key_clears.load(), clear_direct_calls.load());
    swbwa_host_workers_destroy(); swbwa_codec_finish();
    swbwa_mpi_finalize();
    return 0;
}
