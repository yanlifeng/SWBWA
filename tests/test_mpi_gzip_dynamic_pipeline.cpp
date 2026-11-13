#include "swbwa_mpi_gzip.h"
#include "swbwa_pipeline.h"
#include "swbwa_host_workers.h"
#include "swbwa_gzip_output.h"
#include "codec_job.h"
#include "mpi_decoder.h"
extern "C" {
#include "swbwa_mpi.h"
#include "swbwa_output.h"
void kt_pipeline_single(int, void *(*)(void *, int, void *), void *, int);
void kt_pipeline_queue(int, void *(*)(void *, int, void *), void *, int);
}
#include <cassert>
#include <algorithm>
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
#include <sys/uio.h>
#include <unistd.h>

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
        swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
        assert(*static_cast<unsigned *>(item) == f.written++);
        return NULL;
    }
    assert(step == 3);
    swbwa_pipeline_call_on_processor(probe_thread, &f.processor);
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
        assert(f.processor.calls == 5 + 4 * count);
    }
}

static std::string stress;
static int ranks;
static unsigned ready_limit;
static MPI_Comm stress_comm = MPI_COMM_NULL;
static MPI_Request evidence_request = MPI_REQUEST_NULL;
static int evidence = 1;
static bool evidence_sent = false;
static bool alignment_active = false, codec_active = false;
static bool expect_prepared = false;
static thread_local int64_t writing_chunk = -1;
static thread_local uint64_t *prefix_result = NULL;
static thread_local MPI_Win prefix_window = MPI_WIN_NULL;
static thread_local bool flushing_codec = false;
static thread_local bool preparing_output = false;
static bool source_lifecycle_test()
{
    return stress == "source-growth" || stress == "source-deferred";
}
struct DeferredGet {
    void *destination;
    int count, target;
    MPI_Aint displacement;
    MPI_Win window;
};
static thread_local MPI_Win source_window = MPI_WIN_NULL;
static thread_local std::vector<DeferredGet> deferred_gets;
static struct {
    unsigned windows = 0, zero_windows = 0;
    uint64_t gets = 0, completed = 0, flushes = 0, max_pending = 0, max_flush_batch = 0;
} source_observed;
static struct {
    std::atomic<unsigned> pending{0}, overlap_refills{0}, codec_with_pending{0};
    std::atomic<unsigned> decode_with_pending{0}, output_codecs{0}, eof_codecs{0};
    std::atomic<unsigned> prefix_waits{0}, decode_full_queue{0}, prefix_releases{0};
    std::atomic<unsigned> prepared_chunks{0}, prepared_writes{0}, prepare_codecs{0}, freed_sam_chunks{0};
    std::atomic<bool> reader_eof{false};
    unsigned bulk_packets = 0, bulk_payloads = 0, wire_packets = 0, empty_sends = 0, empty_receives = 0;
    unsigned alltoallv_calls = 0, skew_delays = 0;
    std::atomic<uint64_t> decode_calls{0}, decodes_after_overlap{0}, max_overlap_streak{0};
    uint64_t window_decode_calls = 0, overlap_streak = 0;
    std::atomic<uint64_t> read_chunks{0}, processed_chunks{0}, writer_started{0}, written_chunks{0};
    int64_t last_read = -1, last_processed = -1, last_written = -1, last_window = -1;
} observed;

extern "C" ssize_t __real_pwritev(int, const struct iovec *, int, off_t);
extern "C" ssize_t __wrap_pwritev(int fd, const struct iovec *v, int count, off_t offset)
{
    if (stress == "prefix-ahead" && swbwa_mpi_rank() == 0 && writing_chunk == 0 &&
        !observed.prefix_releases.load()) {
        const double deadline = MPI_Wtime() + 30;
        int ready = 0, received = 0;
        // Rank 0's first disk write cannot finish before rank 1 writes ID 3.
        // ID 3 depends on publishing rank 0's ID 2 while ID 0 is still here.
        while (!ready) {
            assert(MPI_Iprobe(1, 83, stress_comm, &ready, MPI_STATUS_IGNORE) == MPI_SUCCESS);
            if (MPI_Wtime() >= deadline) {
                std::fprintf(stderr, "prefix ahead did not progress past blocked pwrite\n");
                MPI_Abort(stress_comm, 1);
            }
            if (!ready) std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
        assert(MPI_Recv(&received, 1, MPI_INT, 1, 83, stress_comm, MPI_STATUS_IGNORE) == MPI_SUCCESS);
        assert(received == 1);
        ++observed.prefix_releases;
    }
    return __real_pwritev(fd, v, count, offset);
}

static uint64_t wire_u64(const unsigned char *&p, const unsigned char *end)
{
    assert(end - p >= 8);
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i) value |= uint64_t(*p++) << (8 * i);
    return value;
}

static void check_wire(const unsigned char *p, size_t size)
{
    const unsigned char *end = p + size;
    uint64_t next_record = 0;
    while (p < end) {
        const uint64_t first = wire_u64(p, end), count = wire_u64(p, end);
        const uint64_t bytes = wire_u64(p, end), spans = wire_u64(p, end);
        assert(first >= next_record && count && count <= UINT64_MAX - first);
        next_record = first + count;
        assert(spans && spans <= uint64_t(end - p) / 16);
        const unsigned char *payload = p + size_t(spans) * 16;
        assert(bytes <= uint64_t(end - payload));
        size_t offset = 0;
        uint64_t lines = 0;
        for (uint64_t j = 0; j < spans; ++j) {
            const uint64_t n = wire_u64(p, end), lfs = wire_u64(p, end);
            assert(n && n <= 4096 && n <= bytes - offset);
            assert(uint64_t(std::count(payload + offset, payload + offset + n, '\n')) == lfs);
            offset += size_t(n); lines += lfs;
        }
        assert(p == payload && offset == bytes && lines / 4 == count && lines % 4 == 0);
        ++observed.wire_packets;
        if (bytes >= (1U << 20)) ++observed.bulk_payloads;
        p = payload + size_t(bytes);
    }
    assert(p == end);
}

extern "C" int __real_MPI_Alltoallv(const void *, const int *, const int *, MPI_Datatype,
                                    void *, const int *, const int *, MPI_Datatype, MPI_Comm);
extern "C" int __wrap_MPI_Alltoallv(const void *send, const int *sc, const int *sd, MPI_Datatype st,
                                    void *receive, const int *rc, const int *rd, MPI_Datatype rt,
                                    MPI_Comm comm)
{
    ++observed.alltoallv_calls;
    assert(st == MPI_BYTE && rt == MPI_BYTE && send && receive);
    size_t sent = 0, received = 0;
    for (int i = 0; i < ranks; ++i) {
        assert(sc[i] >= 0 && rc[i] >= 0 && size_t(sd[i]) == sent && size_t(rd[i]) == received);
        sent += size_t(sc[i]); received += size_t(rc[i]);
        if (sc[i] >= (1 << 20)) ++observed.bulk_packets;
        if (sc[i]) check_wire(static_cast<const unsigned char *>(send) + sd[i], size_t(sc[i]));
    }
    if (!sent) ++observed.empty_sends;
    if (!received) ++observed.empty_receives;
    return __real_MPI_Alltoallv(send, sc, sd, st, receive, rc, rd, rt, comm);
}

extern "C" int __real_MPI_Fetch_and_op(const void *, void *, MPI_Datatype, int,
                                      MPI_Aint, MPI_Op, MPI_Win);
extern "C" int __wrap_MPI_Fetch_and_op(const void *origin, void *result, MPI_Datatype type,
                                      int target, MPI_Aint displacement, MPI_Op op, MPI_Win win)
{
    const int status = __real_MPI_Fetch_and_op(origin, result, type, target, displacement, op, win);
    if (stress == "prefix-full" && writing_chunk == 1 && op == MPI_NO_OP &&
        type == MPI_UINT64_T && target == 0 && displacement == 0) {
        prefix_result = static_cast<uint64_t *>(result);
        prefix_window = win;
    }
    return status;
}

extern "C" int __real_MPI_Win_flush(int, MPI_Win);
extern "C" int __real_MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype, MPI_Win);
extern "C" int __real_MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);
extern "C" int __wrap_MPI_Win_create(void *base, MPI_Aint size, int unit, MPI_Info info,
                                     MPI_Comm comm, MPI_Win *win)
{
    const int status = __real_MPI_Win_create(base, size, unit, info, comm, win);
    // Reader source windows use byte displacements; ticket windows do not.
    if (source_lifecycle_test() && unit == 1) {
        assert(status == MPI_SUCCESS && deferred_gets.empty());
        source_window = *win;
        ++source_observed.windows;
        source_observed.zero_windows += size == 0;
        std::fprintf(stderr, "[test source_window rank=%d] build=%u capacity_bytes=%llu\n",
                     swbwa_mpi_rank(), source_observed.windows, (unsigned long long)size);
    }
    return status;
}

extern "C" int __wrap_MPI_Get(void *destination, int count, MPI_Datatype type, int target,
                              MPI_Aint displacement, int target_count, MPI_Datatype target_type, MPI_Win win)
{
    if (!source_lifecycle_test() || win != source_window)
        return __real_MPI_Get(destination, count, type, target, displacement, target_count, target_type, win);
    assert(type == MPI_BYTE && target_type == MPI_BYTE && count == target_count && count > 0);
    assert(target != swbwa_mpi_rank());
    // Leave both FASTQ and newline metadata unavailable until remote completion.
    std::memset(destination, 0xa5, size_t(count));
    deferred_gets.push_back({destination, count, target, displacement, win});
    ++source_observed.gets;
    source_observed.max_pending = std::max(source_observed.max_pending, uint64_t(deferred_gets.size()));
    return MPI_SUCCESS;
}

extern "C" int __wrap_MPI_Win_flush(int target, MPI_Win win)
{
    uint64_t completed = 0;
    for (const auto &get : deferred_gets) if (get.target == target && get.window == win) {
        const auto *bytes = static_cast<const unsigned char *>(get.destination);
        assert(std::all_of(bytes, bytes + get.count, [](unsigned char c) { return c == 0xa5; }));
        assert(__real_MPI_Get(get.destination, get.count, MPI_BYTE, target, get.displacement,
                              get.count, MPI_BYTE, win) == MPI_SUCCESS);
        ++completed;
    }
    const int status = __real_MPI_Win_flush(target, win);
    if (completed) {
        assert(status == MPI_SUCCESS);
        source_observed.completed += completed;
        ++source_observed.flushes;
        source_observed.max_flush_batch = std::max(source_observed.max_flush_batch, completed);
        deferred_gets.erase(std::remove_if(deferred_gets.begin(), deferred_gets.end(),
            [target, win](const DeferredGet &get) { return get.target == target && get.window == win; }),
            deferred_gets.end());
    }
    if (prefix_result && win == prefix_window && target == 0) {
        if (status == MPI_SUCCESS && *prefix_result < uint64_t(writing_chunk))
            ++observed.prefix_waits;
        prefix_result = NULL;
    }
    return status;
}

extern "C" int __real_swbwa_mpi_decoder_next(swbwa_mpi_decoder *, swbwa_mpi_decoded_batch *);
extern "C" int __wrap_swbwa_mpi_decoder_next(swbwa_mpi_decoder *d, swbwa_mpi_decoded_batch *b)
{
    if (observed.pending.load()) ++observed.decode_with_pending;
    ++observed.decode_calls;
    if (observed.max_overlap_streak >= 2) ++observed.decodes_after_overlap;
    return __real_swbwa_mpi_decoder_next(d, b);
}

extern "C" void __real_swbwa_codec_launch(swbwa_codec_batch *);
static void launch(void *p)
{
    assert(std::this_thread::get_id() == owner);
    assert(!alignment_active && !codec_active);
    codec_active = true;
    swbwa_codec_batch *batch = static_cast<swbwa_codec_batch *>(p);
    const bool output = batch->count && batch->jobs[0].operation == SWBWA_CODEC_GZIP;
    if (output) ++observed.output_codecs;
    else {
        if (observed.pending.load()) ++observed.codec_with_pending;
        if (stress == "prefix-full" && swbwa_mpi_rank() == 1 && observed.prefix_waits.load() &&
            observed.processed_chunks.load() - observed.writer_started.load() == SWBWA_PIPELINE_QUEUE_CAPACITY) {
            ++observed.decode_full_queue;
            if (!evidence_sent) {
                // Release rank 0 only after real prefix polling and decoder service
                // with a full output queue, not after a timing-based delay.
                assert(MPI_Isend(&evidence, 1, MPI_INT, 0, 81, stress_comm, &evidence_request) == MPI_SUCCESS);
                evidence_sent = true;
            }
        }
    }
    __real_swbwa_codec_launch(batch);
    codec_active = false;
}
extern "C" void __wrap_swbwa_codec_launch(swbwa_codec_batch *b)
{
    if (preparing_output) {
        assert(std::this_thread::get_id() == owner);
        assert(b->count && b->jobs[0].operation == SWBWA_CODEC_GZIP);
        ++observed.prepare_codecs;
    }
    swbwa_pipeline_call_on_processor(launch, b);
    if (flushing_codec) ++observed.eof_codecs;
}
extern "C" void __real_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *, unsigned);
extern "C" void __wrap_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *c, unsigned workers)
{
    __real_swbwa_mpi_decoder_config_init(c, workers);
    c->chunk_bytes = 4096; c->jobs_per_rank = 1; c->initial_output_bytes = 65536;
    c->decoded_bytes = 8U << 20; c->memory_bytes = 64U << 20; c->lookahead_bytes = 1U << 20;
    if (stress == "source-pressure") c->decoded_bytes = 65536;
    if (stress == "exchange-bulk" || stress == "exchange-bulk-chunks" || stress == "source-growth") {
        // One compressed member on rank 0 redistributes >1 MiB to each peer.
        c->chunk_bytes = 1U << 20; c->initial_output_bytes = 8U << 20;
    }
}
struct Batch {
    swbwa_mpi_gzip_batch input = {};
    void *prepared = NULL;
    std::vector<std::string> text;
    std::vector<bseq1_t> seqs;
    ~Batch()
    {
        assert(!prepared);
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

static void wait_for_full_queue_decoder()
{
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(15);
    int ready = 0;
    while (!ready) {
        assert(MPI_Iprobe(1, 81, stress_comm, &ready, MPI_STATUS_IGNORE) == MPI_SUCCESS);
        if (std::chrono::steady_clock::now() >= deadline) {
            std::fprintf(stderr, "prefix/full-output decoder gate timed out\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            std::abort();
        }
        if (!ready) std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    int received = 0;
    assert(MPI_Recv(&received, 1, MPI_INT, 1, 81, stress_comm, MPI_STATUS_IGNORE) == MPI_SUCCESS);
    assert(received == 1);
    ++observed.prefix_releases;
}

static void *process(void *context, int step, void *item)
{
    Fixture &f = *static_cast<Fixture *>(context);
    Batch *b = static_cast<Batch *>(item);
    if (step == 0) {
        b = new Batch;
        const int status = swbwa_mpi_gzip_next(f.reader, &b->input);
        assert(deferred_gets.empty());
        assert(status >= 0);
        if (!status) { observed.reader_eof = true; delete b; return NULL; }
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
        assert(!alignment_active && !codec_active);
        alignment_active = true;
        assert(b->input.chunk_id > observed.last_processed);
        observed.last_processed = b->input.chunk_id;
        if ((stress == "source-skew" || source_lifecycle_test()) &&
            (b->input.chunk_id == 0 || f.rank == ranks - 1)) {
            ++observed.skew_delays;
            std::this_thread::sleep_for(std::chrono::milliseconds(8));
        }
        if (!stress.empty() && stress != "source-skew" && !source_lifecycle_test() && f.rank == 0 &&
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
        alignment_active = false;
        if (b->input.records) {
            preparing_output = true;
            assert(!swbwa_output_prepare_chunk(b->seqs.data(), NULL,
                                               int(b->seqs.size()), &b->prepared));
            assert(!swbwa_output_publish_prepared_chunk(b->input.chunk_id,
                                                        int(b->seqs.size()), b->prepared));
            preparing_output = false;
            assert(bool(b->prepared) == expect_prepared);
            if (b->prepared) {
                ++observed.prepared_chunks;
#if !SWBWA_OUTPUT_RMA_ONLY
                // The sink must use compressed storage, not borrow SAM after preparation.
                std::vector<std::string>().swap(b->text);
                for (auto &seq : b->seqs) seq.sam = NULL;
                ++observed.freed_sam_chunks;
#endif
            }
            ++observed.processed_chunks;
        }
        return b;
    }
    if (step == 2) {
        if (b->input.records) {
            assert(b->input.chunk_id > observed.last_written);
            observed.last_written = b->input.chunk_id;
            ++observed.writer_started;
            if (b->input.chunk_id == 0 && stress == "prefix-full") wait_for_full_queue_decoder();
            if (b->input.chunk_id == 0 && (stress == "decoded-refills" || stress == "single-direct"))
                wait_for_overlap();
            if (!stress.empty() && f.rank == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(4));
            if (b->input.chunk_id % 5 == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            writing_chunk = b->input.chunk_id;
            const bseq1_t *seqs = b->seqs.data();
#if !SWBWA_OUTPUT_RMA_ONLY
            if (b->prepared) seqs = NULL;
#endif
            assert(swbwa_output_write_prepared_chunk(b->input.chunk_id, b->input.first_record,
                b->input.first_record + b->input.records, seqs, NULL,
                int(b->seqs.size()), f.mates == 2, b->prepared) == 0);
            if (stress == "prefix-ahead" && f.rank == 1 && b->input.chunk_id == 3) {
                assert(!evidence_sent);
                assert(MPI_Isend(&evidence, 1, MPI_INT, 0, 83, stress_comm, &evidence_request) == MPI_SUCCESS);
                evidence_sent = true;
            }
            if (b->prepared) ++observed.prepared_writes;
            b->prepared = NULL;
            writing_chunk = -1;
            ++observed.written_chunks;
            // Decrement before complete(), which can immediately wake a decoding reader.
            assert(observed.pending.fetch_sub(1) > 0);
            assert(!swbwa_mpi_gzip_complete(f.reader, b->input.chunk_id));
        }
        delete b;
        return NULL;
    }
    assert(step == 3);
    if (stress == "eof-codec") {
        // The writer reaches this callback only after the output queue closes.
        // A real empty-member encode must now be serviced by sync.drain().
        assert(observed.reader_eof.load());
        swbwa_gzip_output *flush = swbwa_gzip_output_create(64);
        assert(flush);
        const struct iovec *vectors = NULL;
        unsigned count = 0;
        uint64_t bytes = 0;
        flushing_codec = true;
        assert(!swbwa_gzip_output_bytes(flush, NULL, 0, &vectors, &count, &bytes));
        flushing_codec = false;
        assert(count == 1 && bytes > 0 && vectors);
        swbwa_gzip_output_destroy(flush);
    }
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
    const char *prepare = std::getenv("SWBWA_GZIP_PREPARE_OUTPUT");
    const std::string output_name(argv[3]);
    expect_prepared = (!prepare || !std::strcmp(prepare, "1")) &&
                      output_name.size() >= 3 && output_name.substr(output_name.size() - 3) == ".gz";
    assert(!swbwa_mpi_init(&argc, &argv));
    ranks = swbwa_mpi_size();
    const char *profile = std::getenv("PIPELINE_TEST_STRESS");
    stress = profile ? profile : "";
    assert(stress.empty() || stress == "decoded-refills" || stress == "window-reuse" ||
           stress == "single-direct" || stress == "window-growth" || stress == "pending-close" ||
           stress == "prefix-full" || stress == "eof-codec" || stress == "drain-mismatch" ||
           stress == "exchange-bulk" || stress == "exchange-bulk-chunks" || stress == "early-first-mismatch" ||
           stress == "prefix-ahead" || stress == "prepare-ahead-mismatch" || stress == "source-skew" ||
           stress == "source-pressure" || source_lifecycle_test());
    if (stress == "drain-mismatch")
        assert(!setenv("SWBWA_GZIP_DECODE_DRAIN", swbwa_mpi_rank() ? "1" : "0", 1));
    if (stress == "early-first-mismatch")
        assert(!setenv("SWBWA_GZIP_EARLY_FIRST", swbwa_mpi_rank() ? "1" : "0", 1));
    if (stress == "prepare-ahead-mismatch")
        assert(!setenv("SWBWA_GZIP_PREPARE_AHEAD", swbwa_mpi_rank() ? "1" : "0", 1));
    if (stress == "prefix-full" || stress == "prefix-ahead") {
        assert(ranks == 2);
        assert(MPI_Comm_dup(MPI_COMM_WORLD, &stress_comm) == MPI_SUCCESS);
    }
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
                                 stress == "window-growth" || stress == "prefix-full" || stress == "prefix-ahead" ||
                                 stress == "exchange-bulk" || stress == "source-growth" ? 1 :
                                 stress == "exchange-bulk-chunks" ? 1500000 :
                                 stress == "source-pressure" ? 65536 :
                                 stress.empty() ? 2048 : 1024, 64, 1);
    assert(f.reader);
    if (stress == "pending-close") reject_pending_close(f.reader, f.rank);
    assert(!swbwa_output_open_online(argv[3], 1));
    assert(!swbwa_mpi_progress_thread_start(0));
    if (std::strcmp(argv[4], "queue") == 0) kt_pipeline_queue(3, process, &f, 3);
    else kt_pipeline_single(1, process, &f, 3);
    assert(observed.pending.load() == 0);
    assert(observed.read_chunks.load() == observed.processed_chunks.load() &&
           observed.read_chunks.load() == observed.written_chunks.load());
    assert(observed.prepared_chunks.load() == (expect_prepared ? observed.written_chunks.load() : 0));
    assert(observed.prepared_chunks.load() == observed.prepared_writes.load() &&
           observed.prepared_chunks.load() == observed.prepare_codecs.load());
#if !SWBWA_OUTPUT_RMA_ONLY
    assert(observed.freed_sam_chunks.load() == observed.prepared_chunks.load());
#endif
    if (stress == "decoded-refills" && !f.rank) {
        assert(observed.max_overlap_streak >= 2);
        assert(observed.decodes_after_overlap > 0);
    }
    if (stress == "single-direct") assert(observed.codec_with_pending.load() > 0);
    if (stress == "prefix-full") {
        if (!f.rank) assert(observed.prefix_releases.load() == 1);
        else {
            assert(observed.prefix_waits.load() && observed.decode_full_queue.load());
            assert(MPI_Wait(&evidence_request, MPI_STATUS_IGNORE) == MPI_SUCCESS);
        }
        assert(MPI_Comm_free(&stress_comm) == MPI_SUCCESS);
    }
    if (stress == "prefix-ahead") {
        if (!f.rank) assert(observed.prefix_releases.load() == 1);
        else {
            assert(evidence_sent);
            assert(MPI_Wait(&evidence_request, MPI_STATUS_IGNORE) == MPI_SUCCESS);
        }
        assert(MPI_Comm_free(&stress_comm) == MPI_SUCCESS);
    }
    if (stress == "eof-codec") assert(observed.eof_codecs.load() == 1);
    if (source_lifecycle_test()) {
        assert(source_observed.gets == source_observed.completed);
        std::fprintf(stderr,
            "[test source_lifecycle rank=%d] windows=%u zero_windows=%u deferred_gets=%llu "
            "completed_gets=%llu deferred_flushes=%llu max_pending_gets=%llu max_flush_batch=%llu\n",
            f.rank, source_observed.windows, source_observed.zero_windows,
            (unsigned long long)source_observed.gets, (unsigned long long)source_observed.completed,
            (unsigned long long)source_observed.flushes, (unsigned long long)source_observed.max_pending,
            (unsigned long long)source_observed.max_flush_batch);
    }
    std::fprintf(stderr,
        "[test pipeline rank=%d] chunks=%llu decode_calls=%llu overlap_refills=%u "
        "max_overlap_streak=%llu decodes_after_overlap=%llu codec_with_pending=%u "
        "decode_with_pending=%u output_codecs=%u eof_codecs=%u prefix_waits=%u "
        "decode_full_queue=%u prefix_releases=%u prepared_chunks=%u prepared_writes=%u "
        "prepare_codecs=%u freed_sam_chunks=%u bulk_packets=%u bulk_payloads=%u wire_packets=%u "
        "empty_sends=%u empty_receives=%u alltoallv_calls=%u skew_delays=%u\n",
        f.rank, (unsigned long long)observed.read_chunks.load(), (unsigned long long)observed.decode_calls.load(),
        observed.overlap_refills.load(), (unsigned long long)observed.max_overlap_streak.load(),
        (unsigned long long)observed.decodes_after_overlap.load(), observed.codec_with_pending.load(),
        observed.decode_with_pending.load(), observed.output_codecs.load(), observed.eof_codecs.load(),
        observed.prefix_waits.load(), observed.decode_full_queue.load(), observed.prefix_releases.load(),
        observed.prepared_chunks.load(), observed.prepared_writes.load(),
        observed.prepare_codecs.load(), observed.freed_sam_chunks.load(),
        observed.bulk_packets, observed.bulk_payloads, observed.wire_packets,
        observed.empty_sends, observed.empty_receives, observed.alltoallv_calls, observed.skew_delays);
    assert(!swbwa_output_ordered_finish(swbwa_mpi_gzip_chunk_count(f.reader)));
    assert(!swbwa_mpi_gzip_report(f.reader, stderr));
    swbwa_mpi_gzip_close(f.reader);
    assert(!swbwa_output_close());
    assert(!swbwa_mpi_progress_thread_stop());
    check_handoff_reuse();
    const unsigned pipelines = 2 + unsigned(std::strcmp(argv[4], "queue") == 0);
    const char *ahead = std::getenv("SWBWA_GZIP_PREPARE_AHEAD");
    const unsigned producers = ahead && !std::strcmp(ahead, "1") && ranks > 1 &&
        SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC && !std::strcmp(argv[4], "queue");
    assert(key_sets.load() == 2 * pipelines + producers &&
           key_clears.load() == 2 * pipelines + producers &&
           clear_direct_calls.load() == 2 * pipelines + producers);
    std::fprintf(stderr,
        "[test pipeline keys rank=%d] pipelines=%u sets=%u clears=%u clear_direct_calls=%u\n",
        f.rank, pipelines, key_sets.load(), key_clears.load(), clear_direct_calls.load());
    swbwa_host_workers_destroy(); swbwa_codec_finish();
    swbwa_mpi_finalize();
    return 0;
}
