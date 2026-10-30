// Incremental form of SWBWA's distributed gzip protocol experiment. Only
// metadata, CRCs and bounded dictionaries cross ranks; decoded bodies do not.
#include "mpi_decoder.h"
#include "mpi_decoder_protocol.h"
#include "codec_job.h"
#include "crc_combine.h"
#include "parallel_pread.h"
#include <zlib.h>

#include <cerrno>
#include <atomic>
#include <chrono>
#include <climits>
#include <condition_variable>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
#include <athread.h>
#endif

namespace {
using namespace swbwa_mpi_gzip;
constexpr size_t GUARD = 128U << 10;
constexpr size_t HEADER_LIMIT = 1U << 20;

void mpi_check(int status) {
    if (status == MPI_SUCCESS) return;
    char text[MPI_MAX_ERROR_STRING];
    int length = 0;
    MPI_Error_string(status, text, &length);
    throw std::runtime_error(std::string(text, length));
}

std::runtime_error io_error(const char *operation) {
    return std::runtime_error(std::string(operation) + ": " + std::strerror(errno));
}

class Source {
public:
    Source(const char *path, swbwa_mpi_decoder_metrics &metrics) : metrics_(metrics) {
        fd_ = ::open(path, O_RDONLY);
        if (fd_ < 0) throw io_error("open gzip");
        try {
            struct stat info{};
            if (fstat(fd_, &info)) throw io_error("stat gzip");
            if (!S_ISREG(info.st_mode) || info.st_size < 18 ||
                uint64_t(info.st_size) > UINT64_MAX / 8)
                throw std::runtime_error("input must be a regular gzip file of at least 18 bytes");
            size_ = uint64_t(info.st_size);
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
            static std::atomic<unsigned> lane{0};
            unsigned long home;
            __asm__ volatile("rcid %0" : "=r"(home));
            io_target_ = int(((home & 7) + 1 + lane.fetch_add(1) %
                (SWBWA_HOST_MPE_THREADS - 1)) % SWBWA_HOST_MPE_THREADS);
#endif
            metrics_.prefetch_cg = io_target_;
            unsigned readers = 1;
#if SWBWA_HOST_MPE_THREADS > 1
            readers = SWBWA_HOST_MPE_THREADS;
#endif
            if (const char *setting = std::getenv("SWBWA_GZIP_INPUT_READERS")) {
                char *end = nullptr;
                errno = 0;
                const unsigned long value = std::strtoul(setting, &end, 10);
                if (errno || !*setting || *end || value < 1 || value > readers)
                    throw std::runtime_error("gzip input readers exceeds available MPEs");
                readers = unsigned(value);
            }
            reads_.reset(new SwbwaParallelPread(readers));
            if (const char *setting = std::getenv("SWBWA_GZIP_PREFETCH_MIB")) {
                char *end = nullptr;
                errno = 0;
                unsigned long value = std::strtoul(setting, &end, 10);
                if (errno || !*setting || *end || value > 256)
                    throw std::runtime_error("gzip prefetch must be 0..256 MiB");
                prefetch_limit_ = size_t(value) << 20;
            }
        } catch (...) { ::close(fd_); fd_ = -1; throw; }
    }
    ~Source() {
        if (pending_.valid()) pending_.wait();
        if (io_worker_.joinable()) {
            {
                std::lock_guard<std::mutex> lock(io_mutex_);
                io_stopping_ = true;
            }
            io_ready_.notify_one();
            io_worker_.join();
        }
        if (fd_ >= 0) ::close(fd_);
    }
    uint64_t size() const { return size_; }
    void read(uint64_t offset, void *buffer, size_t count) {
        if (offset > size() || count > size() - offset)
            throw std::runtime_error("truncated input");
        const double begin = MPI_Wtime();
        Io io;
        read_exact(fd_, offset, buffer, count, io);
        metrics_.source_read_bytes += count;
        metrics_.source_read_calls += io.calls;
        metrics_.source_io_s += io.seconds;
        metrics_.io_s += MPI_Wtime() - begin;
    }
    std::unique_ptr<unsigned char[]> read_batch(uint64_t offset, size_t count) {
        const double begin = MPI_Wtime();
        Transfer transfer;
        if (pending_.valid()) {
            if (offset != pending_offset_ || count != pending_size_)
                throw std::logic_error("gzip prefetch range disagrees with decoder plan");
            transfer = pending_.get();
            ++metrics_.prefetch_batches;
            metrics_.prefetch_bytes += transfer.loaded - transfer.cached;
            metrics_.prefetch_io_s += transfer.io.seconds;
        } else {
            transfer.data.reset(new unsigned char[count + GUARD]);
        }
        bulk_read(offset + transfer.loaded, transfer.data.get() + transfer.loaded,
                  count - transfer.loaded, transfer.io);
        std::memset(transfer.data.get() + count, 0, GUARD);
        metrics_.source_read_bytes += count - transfer.cached;
        metrics_.source_cache_bytes += transfer.cached;
        metrics_.source_read_calls += transfer.io.calls;
        metrics_.source_io_s += transfer.io.seconds;
        metrics_.io_s += MPI_Wtime() - begin;
        return std::move(transfer.data);
    }
    void prefetch(uint64_t offset, size_t count, const unsigned char *previous = nullptr,
                  uint64_t previous_offset = 0, size_t previous_size = 0) {
        if (!prefetch_limit_ || !count) return;
        if (pending_.valid()) throw std::logic_error("unconsumed gzip prefetch request");
        pending_offset_ = offset; pending_size_ = count;
        const size_t loaded = std::min(count, prefetch_limit_);
        auto transfer = std::make_shared<Transfer>();
        transfer->data.reset(new unsigned char[count + GUARD]);
        const uint64_t shared_begin = std::max(offset, previous_offset);
        const uint64_t shared_end = std::min(offset + loaded, previous_offset + previous_size);
        size_t prefix = loaded, suffix = loaded;
        const char *reuse = std::getenv("SWBWA_GZIP_REUSE_INPUT");
        if (previous && (!reuse || std::strcmp(reuse, "0")) && shared_begin < shared_end) {
            prefix = size_t(shared_begin - offset);
            suffix = size_t(shared_end - offset);
            transfer->cached = suffix - prefix;
            const auto begin = std::chrono::steady_clock::now();
            std::memcpy(transfer->data.get() + prefix, previous + (shared_begin - previous_offset),
                        transfer->cached);
            metrics_.source_cache_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        }
        // Only the exact next rank-local stripe is read. This worker never
        // calls MPI or enters the CPE runtime, and hands its allocation over.
        std::packaged_task<Transfer()> request([this, offset, loaded, prefix, suffix, transfer]() {
            // Bind inside the task so binding failures also reach the future.
            if (!io_bound_) {
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
                SET_TO_SWCG(io_target_);
                unsigned long actual;
                __asm__ volatile("rcid %0" : "=r"(actual));
                if (int(actual & 7) != io_target_)
                    throw std::runtime_error("gzip I/O worker MPE binding failed");
#endif
                io_bound_ = true;
            }
            bulk_read(offset, transfer->data.get(), prefix, transfer->io);
            bulk_read(offset + suffix, transfer->data.get() + suffix, loaded - suffix, transfer->io);
            transfer->loaded = loaded;
            return std::move(*transfer);
        });
        auto pending = request.get_future();
        {
            std::lock_guard<std::mutex> lock(io_mutex_);
            if (!io_worker_.joinable())
                io_worker_ = std::thread(&Source::prefetch_loop, this);
            io_request_ = std::move(request);
            pending_ = std::move(pending);
        }
        io_ready_.notify_one();
    }
    uint64_t header(uint64_t begin, uint32_t &fingerprint) {
        // Buffer only the header cursor, never scan the compressed body.
        std::array<unsigned char, 4096> buffer{};
        uint64_t pos = begin, loaded = begin;
        size_t available = 0;
        uint32_t crc = 0;
        auto take = [&]() {
            if (pos - begin >= HEADER_LIMIT) throw std::runtime_error("gzip header exceeds 1 MiB");
            if (pos == loaded + available) {
                loaded = pos;
                available = size_t(std::min<uint64_t>(buffer.size(), size() - pos));
                if (!available) throw std::runtime_error("truncated gzip header");
                read(pos, buffer.data(), available);
            }
            unsigned char c = buffer[size_t(pos++ - loaded)];
            crc = crc32(crc, &c, 1);
            return unsigned(c);
        };
        if (take() != 31 || take() != 139 || take() != 8)
            throw std::runtime_error("invalid gzip header/method");
        unsigned flags = take();
        if (flags & 0xe0) throw std::runtime_error("reserved gzip header flags");
        for (unsigned i = 0; i < 6; ++i) take();
        if (flags & 4) {
            unsigned n = take(); n |= take() << 8;
            for (unsigned i = 0; i < n; ++i) take();
        }
        for (unsigned flag : {8U, 16U}) if (flags & flag) while (take()) {}
        fingerprint = crc;
        if (flags & 2) {
            unsigned stored = take(); stored |= take() << 8;
            if (stored != (fingerprint & 65535)) throw std::runtime_error("gzip header CRC mismatch");
        }
        if (pos > size() - 8) throw std::runtime_error("truncated gzip header");
        return pos * 8;
    }
    uint64_t trailer_word(uint64_t offset) {
        unsigned char bytes[8];
        read(offset, bytes, sizeof(bytes));
        uint64_t word = 0;
        for (unsigned i = 0; i < 8; ++i) word |= uint64_t(bytes[i]) << (8 * i);
        return word;
    }
    uint64_t validate_member(uint64_t cursor, uint32_t crc, uint64_t bytes) {
        const uint64_t offset = (cursor + 7) / 8;
        const uint64_t end = add(offset, 8);
        if (end > size()) throw std::runtime_error("truncated trailer");
        uint64_t trailer = trailer_word(offset);
        if (uint32_t(trailer) != crc || uint32_t(trailer >> 32) != uint32_t(bytes))
            throw std::runtime_error("gzip member CRC/ISIZE mismatch");
        return end;
    }
private:
    struct Io { uint64_t calls = 0; double seconds = 0; };
    struct Transfer {
        std::unique_ptr<unsigned char[]> data;
        size_t loaded = 0;
        size_t cached = 0;
        Io io;
    };
    void bulk_read(uint64_t offset, void *buffer, size_t count, Io &io) {
        const auto result = reads_->read(fd_, offset, buffer, count);
        io.calls += result.calls;
        io.seconds += result.seconds;
    }
    void prefetch_loop() {
        for (;;) {
            std::packaged_task<Transfer()> request;
            {
                std::unique_lock<std::mutex> lock(io_mutex_);
                io_ready_.wait(lock, [this] { return io_stopping_ || io_request_.valid(); });
                if (io_stopping_) return;
                request = std::move(io_request_);
            }
            request();
        }
    }
    static void read_exact(int fd, uint64_t offset, void *buffer, size_t count, Io &io) {
        auto begin = std::chrono::steady_clock::now();
        unsigned char *p = static_cast<unsigned char *>(buffer);
        while (count) {
            ++io.calls;
            ssize_t n = pread(fd, p, count, static_cast<off_t>(offset));
            if (n < 0 && errno == EINTR) continue;
            if (n < 0) throw io_error("pread gzip");
            if (!n) throw std::runtime_error("truncated or changed input");
            count -= size_t(n); p += n; offset += uint64_t(n);
        }
        io.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    }
    swbwa_mpi_decoder_metrics &metrics_;
    int fd_ = -1;
    uint64_t size_ = 0;
    size_t prefetch_limit_ = size_t(128) << 20;
    uint64_t pending_offset_ = 0;
    size_t pending_size_ = 0;
    std::future<Transfer> pending_;
    int io_target_ = -1;
    std::mutex io_mutex_;
    std::condition_variable io_ready_;
    std::packaged_task<Transfer()> io_request_;
    std::thread io_worker_;
    bool io_stopping_ = false;
    bool io_bound_ = false; // Accessed only by the persistent worker.
    std::unique_ptr<SwbwaParallelPread> reads_;
};

template<class T> void release(std::vector<T> &v) { std::vector<T>().swap(v); }

class AccumulateTime {
public:
    explicit AccumulateTime(double &seconds) : seconds_(seconds), begin_(MPI_Wtime()) {}
    ~AccumulateTime() { seconds_ += MPI_Wtime() - begin_; }
private:
    double &seconds_;
    const double begin_;
};
} // namespace

struct swbwa_mpi_decoder {
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0, ranks = 0;
    swbwa_mpi_decoder_config config{};
    swbwa_mpi_decoder_metrics metrics{};
    std::unique_ptr<Source> source;
    char error[512]{};
    uint64_t cursor = 0, decoded = 0, epoch = 0;
    uint32_t crc = 0, member_crc = 0;
    uint64_t member_begin = 0;
    bool eof = false, failed = false;
    size_t capacity = 0;
    uint64_t fixed_bytes = 0;
    Window history;
    SwbwaCrcCombine combine;
    std::unique_ptr<unsigned char[]> input;
    size_t input_size = 0;
    uint64_t input_base = 0;
    unsigned gap_repair_limit = 8;
    bool parallel_history = false;
    std::unique_ptr<unsigned char[]> bytes;
    size_t bytes_size = 0;
    std::vector<std::unique_ptr<uint16_t[]>> symbols;
    std::vector<std::unique_ptr<unsigned char[]>> tails;
    std::vector<std::vector<uint32_t>> lines;
    std::vector<swbwa_mpi_decoded_extent> extents;

    void init(MPI_Comm parent, const char *path, const swbwa_mpi_decoder_config *options) {
        AccumulateTime timer(metrics.open_s);
        int inter = 0;
        mpi_check(MPI_Comm_test_inter(parent, &inter));
        if (inter) throw std::runtime_error("decoder requires an intracommunicator");
        mpi_check(MPI_Comm_dup(parent, &comm));
        mpi_check(MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN));
        mpi_check(MPI_Comm_rank(comm, &rank));
        mpi_check(MPI_Comm_size(comm, &ranks));
        if (!path || !*path || std::strlen(path) > 4096)
            throw std::runtime_error("invalid gzip path");
        if (options) config = *options;
        else swbwa_mpi_decoder_config_init(&config, 384);
        if (const char *setting = std::getenv("SWBWA_GZIP_GAP_REPAIR_LIMIT")) {
            char *end = nullptr;
            errno = 0;
            const unsigned long value = std::strtoul(setting, &end, 10);
            if (errno || !*setting || *end || value > 8)
                throw std::runtime_error("gzip gap repair limit must be 0..8");
            gap_repair_limit = unsigned(value);
        }
        if (const char *setting = std::getenv("SWBWA_GZIP_PARALLEL_HISTORY")) {
            if (std::strcmp(setting, "0") && std::strcmp(setting, "1"))
                throw std::runtime_error("SWBWA_GZIP_PARALLEL_HISTORY must be 0 or 1");
            parallel_history = setting[0] == '1';
        }
        if (config.chunk_bytes < 1024 || config.chunk_bytes > (64U << 20) ||
            config.lookahead_bytes < 1024 || config.lookahead_bytes > (64U << 20) ||
            config.memory_bytes < (2U << 20) || config.decoded_bytes < 65536 ||
            config.initial_output_bytes < 65536 || config.initial_output_bytes > UINT_MAX ||
            !config.jobs_per_rank || config.jobs_per_rank > 384 ||
            !config.workers || config.workers > 384 || ranks <= 0)
            throw std::runtime_error("invalid decoder configuration");
        const uint64_t slots = multiply(unsigned(ranks), config.jobs_per_rank);
        if (slots > INT_MAX / FIELDS)
            throw std::runtime_error("MPI job metadata exceeds int count limits");
        capacity = config.initial_output_bytes;
        // Conservative adapter accounting: all simultaneous wire arrays, job
        // descriptors and window temporaries, then per-job payload allocations.
        fixed_bytes = add(multiply(slots, sizeof(Meta) + sizeof(Plan) + 32),
            add(multiply(config.jobs_per_rank, sizeof(swbwa_codec_job) * 3 + 512), HISTORY * 4));
        if (parallel_history) fixed_bytes = add(fixed_bytes,
            add(multiply(unsigned(ranks), sizeof(HistorySummary)),
                add(HISTORY * 4, HISTORY_TRACE_SLOTS * sizeof(HistoryTraceEntry))));
        metrics.initial_batch_jobs = job_count(capacity);
        if (!metrics.initial_batch_jobs) throw std::runtime_error("initial job does not fit memory/decoded budget");

        std::array<char, 4097> root_path{};
        if (rank == 0) std::memcpy(root_path.data(), path, std::strlen(path));
        mpi_check(MPI_Bcast(root_path.data(), root_path.size(), MPI_CHAR, 0, comm));
        if (std::strcmp(path, root_path.data())) throw std::runtime_error("ranks disagree on input path");
        source.reset(new Source(path, metrics));
        uint32_t header_crc = 0;
        cursor = source->header(0, header_crc);
        std::array<uint64_t, 14> values{{2, source->size(), cursor, header_crc, source->trailer_word(source->size() - 8),
            config.chunk_bytes, config.lookahead_bytes,
            config.memory_bytes, config.decoded_bytes, config.initial_output_bytes,
            config.jobs_per_rank, config.workers, gap_repair_limit, uint64_t(parallel_history)}};
        decltype(values) low{}, high{};
        mpi_check(MPI_Allreduce(values.data(), low.data(), values.size(), MPI_UINT64_T, MPI_MIN, comm));
        mpi_check(MPI_Allreduce(values.data(), high.data(), values.size(), MPI_UINT64_T, MPI_MAX, comm));
        if (low != high) {
            const char *names[] = {"protocol", "file_size", "header_end", "header_crc",
                "trailer", "chunk_bytes",
                "lookahead_bytes", "memory_bytes", "decoded_bytes", "initial_output_bytes",
                "jobs_per_rank", "workers", "gap_repair_limit", "parallel_history"};
            for (size_t i = 0; i < values.size(); ++i) {
                if (low[i] == high[i]) continue;
                char detail[256];
                std::snprintf(detail, sizeof(detail),
                    "ranks disagree on input metadata/configuration: %s local=%llu min=%llu max=%llu",
                    names[i], (unsigned long long)values[i],
                    (unsigned long long)low[i], (unsigned long long)high[i]);
                throw std::runtime_error(detail);
            }
        }
        metrics.open_io_s = metrics.io_s;
    }

    unsigned job_count(size_t cap) const {
        // A next-window prefetch can coexist with the current compressed input.
        const uint64_t overhead = add(fixed_bytes, multiply(2, add(config.lookahead_bytes, GUARD)));
        if (overhead >= config.memory_bytes || cap > config.decoded_bytes) return 0;
        const uint64_t per_job = add(multiply(cap, 3), add(multiply(config.chunk_bytes, 2),
            add(HISTORY * 2, multiply((cap + PAGE - 1) / PAGE, sizeof(uint32_t)))));
        return unsigned(std::min<uint64_t>(config.jobs_per_rank,
            std::min<uint64_t>((config.memory_bytes - overhead) / per_job, config.decoded_bytes / cap)));
    }

    struct InputWindow { unsigned count; uint64_t base; size_t bytes; };
    InputWindow input_window() const {
        unsigned count = job_count(capacity);
        if (!count) throw std::runtime_error("anchored DISCOVER exhausted memory/decoded budget");
        const uint64_t remaining = source->size() - cursor / 8;
        const uint64_t strides = remaining / config.chunk_bytes + (remaining % config.chunk_bytes != 0);
        const uint64_t per_rank = strides / unsigned(ranks) + (strides % unsigned(ranks) != 0);
        count = unsigned(std::min<uint64_t>(count, per_rank));
        if (!count) throw std::runtime_error("no compressed data before validated EOF");
        const unsigned anchor = unsigned(epoch % unsigned(ranks));
        const uint64_t group = (unsigned(rank) + unsigned(ranks) - anchor) % unsigned(ranks);
        const uint64_t base = add(cursor / 8, multiply(multiply(group, count), config.chunk_bytes));
        const size_t bytes = base < source->size() ? size_t(std::min<uint64_t>(source->size() - base,
            add(multiply(count, config.chunk_bytes), config.lookahead_bytes))) : 0;
        return {count, base, bytes};
    }

    void launch(std::vector<swbwa_codec_job> &jobs, bool resolve) {
        if (jobs.empty()) return;
        uint64_t &peak = resolve ? metrics.max_resolve_jobs : metrics.max_discover_jobs;
        peak = std::max<uint64_t>(peak, jobs.size());
        swbwa_codec_batch batch{};
        batch.jobs = jobs.data(); batch.count = unsigned(jobs.size()); batch.workers = config.workers;
        AccumulateTime timer(resolve ? metrics.resolve_wall_s : metrics.discover_wall_s);
        swbwa_codec_launch(&batch);
        if (resolve) {
            metrics.resolve_run_s += batch.run_s; metrics.resolve_wait_s += batch.wait_s;
        } else {
            metrics.discover_run_s += batch.run_s; metrics.discover_wait_s += batch.wait_s;
            metrics.discover_jobs += jobs.size();
        }
    }

    void describe(const swbwa_codec_job &job, Meta &m) {
        m[OK] = job.status == Z_OK; m[OUTPUT_FULL] = job.output_full != 0;
        m[KNOWN] = job.known_start;
        metrics.discover_failed_jobs += !m[OK];
        metrics.discover_output_full_jobs += m[OUTPUT_FULL];
        if (m[OK]) {
            m[START] = add(m[BASE], job.start_bit);
            m[START_MAX] = add(m[BASE], job.start_max_bit);
            m[END] = add(m[BASE], job.end_bit);
            m[BYTES] = job.output_size; m[FINAL] = job.final_block != 0;
        }
    }

    void repair_gaps(std::vector<swbwa_codec_job> &jobs, std::vector<Meta> &local_meta,
                     std::vector<Meta> &all_meta, std::vector<Plan> &plan,
                     Chain &chain, unsigned count, unsigned anchor) {
        if (!gap_repair_limit) return;
        AccumulateTime timer(metrics.gap_repair_s);
        for (unsigned attempt = 0; attempt < gap_repair_limit; ++attempt) {
            std::array<uint64_t, 2> request{};
            double phase = MPI_Wtime();
            if (!rank) {
                const size_t id = repair_slot(all_meta, plan, chain, unsigned(ranks), count,
                    anchor, cursor / 8 * 8, config.chunk_bytes * 8);
                if (id != all_meta.size()) request = {{id + 1, chain.cursor}};
            }
            metrics.metadata_local_s += MPI_Wtime() - phase;
            phase = MPI_Wtime();
            mpi_check(MPI_Bcast(request.data(), request.size(), MPI_UINT64_T, 0, comm));
            metrics.metadata_mpi_s += MPI_Wtime() - phase;
            if (!request[0]) break;
            const size_t id = size_t(request[0] - 1), j = id % count;
            const int owner = int(id / count);
            Meta replacement{};
            if (rank == owner) {
                phase = MPI_Wtime();
                ++metrics.discarded_jobs;
                if (local_meta[j][OK]) metrics.discarded_output_bytes += local_meta[j][BYTES];
                std::vector<swbwa_codec_job> retry(1, jobs[j]);
                retry[0].start_bit = size_t(request[1] - local_meta[j][BASE]);
                retry[0].known_start = 2;
                retry[0].dictionary = nullptr; retry[0].dictionary_size = 0;
                metrics.prepare_s += MPI_Wtime() - phase;
                // The prefix proves this bit position, but its dictionary is
                // still unknown. Decode to markers; normal RESOLVE supplies it.
                launch(retry, false);
                phase = MPI_Wtime();
                jobs[j] = retry[0];
                describe(jobs[j], local_meta[j]);
                replacement = local_meta[j];
                ++metrics.gap_repair_jobs;
                metrics.gap_repair_ok += replacement[OK];
                metrics.metadata_local_s += MPI_Wtime() - phase;
            }
            phase = MPI_Wtime();
            mpi_check(MPI_Bcast(replacement.data(), FIELDS, MPI_UINT64_T, owner, comm));
            metrics.metadata_mpi_s += MPI_Wtime() - phase;
            phase = MPI_Wtime();
            if (!rank) {
                all_meta[id] = replacement;
                chain = select(all_meta, plan, unsigned(ranks), count, anchor,
                               cursor, decoded, source->size() * 8, capacity);
            }
            metrics.metadata_local_s += MPI_Wtime() - phase;
            // A very long block can require more capacity. Commit the proven
            // prefix and let the regular next-epoch anchor grow its budget.
            if (!replacement[OK]) break;
        }
    }

    void clear_payload() {
        AccumulateTime timer(metrics.cleanup_s);
        release(extents); bytes.reset(); bytes_size = 0;
        release(lines); release(symbols); release(tails); input.reset(); input_size = 0;
    }

    int next(swbwa_mpi_decoded_batch &out) {
        AccumulateTime timer(metrics.next_s);
        if (failed) throw std::runtime_error("decoder failed; caller must MPI_Abort");
        clear_payload();
        out.epoch = epoch; out.global_begin = decoded;
        out.global_decoded_end = decoded; out.stream_crc32 = crc; out.eof = eof;
        if (eof) {
            // Remain collective even on repeated EOF calls.
            AccumulateTime eof_timer(metrics.eof_barrier_s);
            mpi_check(MPI_Barrier(comm));
            return 0;
        }
        double phase = MPI_Wtime();
        const unsigned anchor = unsigned(epoch % unsigned(ranks));
        unsigned count = 0;
        std::vector<swbwa_codec_job> jobs;
        std::vector<Meta> local_meta, all_meta;
        metrics.prepare_s += MPI_Wtime() - phase;
        for (;;) {
            const double attempt_begin = MPI_Wtime();
            // A retry must release the entire old symbol/input pool first.
            phase = MPI_Wtime();
            release(jobs); release(local_meta); release(all_meta);
            metrics.cleanup_s += MPI_Wtime() - phase;
            clear_payload();
            phase = MPI_Wtime();
            const double io_before = metrics.io_s;
            const InputWindow window = input_window();
            count = window.count;
            const uint64_t base = window.base;
            input_base = base;
            local_meta.resize(count);
            jobs.reserve(count); symbols.reserve(count); tails.reserve(count);
            if (rank == 0) all_meta.resize(size_t(ranks) * count);
            if (base < source->size()) {
                const size_t n = window.bytes;
                input_size = n + GUARD;
                input = source->read_batch(base, n);
                for (unsigned j = 0; j < count; ++j) {
                    const uint64_t offset = multiply(j, config.chunk_bytes);
                    if (offset >= source->size() - base) break;
                    symbols.emplace_back(new uint16_t[capacity]);
                    tails.emplace_back(new unsigned char[HISTORY]);
                    swbwa_codec_job job{};
                    job.input = input.get() + offset;
                    job.logical_size = size_t(std::min<uint64_t>(source->size() - base - offset,
                        config.chunk_bytes + config.lookahead_bytes));
                    job.input_size = job.logical_size + GUARD;
                    job.output_capacity = capacity;
                    job.symbols = symbols.back().get();
                    job.tail = tails.back().get();
                    job.known_start = unsigned(rank) == anchor && j == 0;
                    job.start_bit = job.known_start ? cursor % 8 : 0;
                    job.stop_bit = config.chunk_bytes * 8;
                    if (job.known_start) {
                        job.dictionary = history.data(); job.dictionary_size = unsigned(history.size());
                    }
                    job.operation = SWBWA_CODEC_DISCOVER;
                    jobs.push_back(job);
                    local_meta[j][PRESENT] = 1;
                    local_meta[j][BASE] = multiply(add(base, offset), 8);
                    local_meta[j][LIMIT] = add(local_meta[j][BASE], job.logical_size * 8);
                    local_meta[j][KNOWN] = job.known_start;
                }
            }
            const uint64_t estimate = add(fixed_bytes, add(input_size, multiply(jobs.size(),
                add(multiply(capacity, 3), add(HISTORY * 2, ((capacity + PAGE - 1) / PAGE) * 4)))));
            if (estimate > config.memory_bytes) throw std::runtime_error("adapter buffer budget exceeded");
            metrics.peak_buffer_bytes = std::max(metrics.peak_buffer_bytes, estimate);
            metrics.prepare_s += MPI_Wtime() - phase - (metrics.io_s - io_before);
            launch(jobs, false);
            phase = MPI_Wtime();
            for (size_t j = 0; j < jobs.size(); ++j) {
                describe(jobs[j], local_meta[j]);
            }
            metrics.metadata_local_s += MPI_Wtime() - phase;
            phase = MPI_Wtime();
            mpi_check(MPI_Gather(local_meta.data(), count * FIELDS, MPI_UINT64_T,
                all_meta.data(), count * FIELDS, MPI_UINT64_T, 0, comm));
            metrics.metadata_mpi_s += MPI_Wtime() - phase;
            phase = MPI_Wtime();
            std::array<uint64_t, 3> retry{};
            if (rank == 0 && all_meta[size_t(anchor) * count][OUTPUT_FULL]) {
                size_t low = capacity, high = size_t(std::min<uint64_t>(UINT_MAX, config.decoded_bytes));
                // Find the largest capacity fitting at least one job, including
                // non-power-of-two budgets, then grow towards it geometrically.
                while (low < high) {
                    size_t middle = low + (high - low + 1) / 2;
                    if (job_count(middle)) low = middle;
                    else high = middle - 1;
                }
                const size_t enlarged = size_t(std::min<uint64_t>(low, multiply(capacity, 2)));
                if (enlarged <= capacity)
                    throw std::runtime_error("anchored DISCOVER exhausted memory/decoded budget");
                retry = {1, enlarged, job_count(enlarged)};
            }
            metrics.metadata_local_s += MPI_Wtime() - phase;
            phase = MPI_Wtime();
            mpi_check(MPI_Bcast(retry.data(), retry.size(), MPI_UINT64_T, 0, comm));
            metrics.metadata_mpi_s += MPI_Wtime() - phase;
            if (!retry[0]) break;
            capacity = size_t(retry[1]);
            if (job_count(capacity) != retry[2]) throw std::runtime_error("retry job grid mismatch");
            ++metrics.capacity_retries;
            metrics.retry_jobs += jobs.size();
            metrics.retry_s += MPI_Wtime() - attempt_begin;
        }

        phase = MPI_Wtime();
        std::vector<Plan> plan(size_t(ranks) * count);
        Chain chain{{}, cursor, decoded, false};
        if (rank == 0) chain = select(all_meta, plan, unsigned(ranks), count, anchor,
                                     cursor, decoded, source->size() * 8, capacity);
        metrics.metadata_local_s += MPI_Wtime() - phase;
        repair_gaps(jobs, local_meta, all_meta, plan, chain, count, anchor);
        phase = MPI_Wtime();
        mpi_check(MPI_Bcast(plan.data(), plan.size() * 3, MPI_UINT64_T, 0, comm));
        metrics.metadata_mpi_s += MPI_Wtime() - phase;

        phase = MPI_Wtime();
        std::vector<std::array<unsigned char, HISTORY>> dictionaries(jobs.size());
        std::vector<swbwa_codec_job> resolved;
        std::vector<size_t> selected;
        resolved.reserve(jobs.size()); selected.reserve(jobs.size()); extents.reserve(jobs.size());
        Window next_history = history;
        size_t local_bytes = 0;
        metrics.prepare_s += MPI_Wtime() - phase;
        const double window_begin = MPI_Wtime(), window_mpi_before = metrics.window_mpi_s;
        const auto prepare_local = [&] {
            for (size_t j = 0; j < jobs.size(); ++j) {
                if (!plan[size_t(rank) * count + j][0]) continue;
                auto job = jobs[j];
                std::copy(next_history.begin(), next_history.end(), dictionaries[j].begin());
                job.dictionary = dictionaries[j].data(); job.dictionary_size = unsigned(next_history.size());
                if (job.tail_known) {
                    next_history.assign(job.tail, job.tail + HISTORY);
                    ++metrics.known_tail_jobs;
                } else {
                    next_history = advance(next_history, symbols[j].get(), job.output_size);
                    ++metrics.repaired_tail_jobs;
                }
                job.operation = SWBWA_CODEC_RESOLVE;
                resolved.push_back(job); selected.push_back(j);
                local_bytes = size_t(add(local_bytes, job.output_size));
            }
        };
        if (parallel_history) {
            std::vector<HistorySummary> summaries(ranks);
            auto &own = summaries[rank];
            std::vector<SymbolSpan> spans;
            spans.reserve(jobs.size());
            for (size_t j = 0; j < jobs.size(); ++j) {
                if (!plan[size_t(rank) * count + j][0]) continue;
                spans.push_back({symbols[j].get(), jobs[j].output_size});
            }
            summarize_history(spans, own);
            // One symbolic transfer function per rank replaces the serialized
            // broadcasts of every rank's locally reconstructed dictionaries.
            phase = MPI_Wtime();
            mpi_check(MPI_Allgather(MPI_IN_PLACE, 0, MPI_BYTE, summaries.data(),
                                    sizeof(HistorySummary), MPI_BYTE, comm));
            metrics.window_mpi_s += MPI_Wtime() - phase;
            Window global_history = history, incoming;
            for (unsigned group = 0; group < unsigned(ranks); ++group) {
                const unsigned owner = (anchor + group) % unsigned(ranks);
                if (summaries[owner].bytes > config.decoded_bytes)
                    throw std::runtime_error("history summary exceeds decoded budget");
                if (unsigned(rank) == owner) incoming = global_history;
                global_history = resolve_history(global_history, summaries[owner]);
            }
            next_history.swap(incoming);
            prepare_local();
            next_history.swap(global_history);
            ++metrics.parallel_history_epochs;
        } else {
            // General case: rank blocks are contiguous in stream order.
            for (unsigned group = 0; group < unsigned(ranks); ++group) {
                const unsigned owner = (anchor + group) % unsigned(ranks);
                bool active = false;
                for (unsigned j = 0; j < count; ++j) active |= plan[size_t(owner) * count + j][0] != 0;
                if (!active) continue;
                if (unsigned(rank) == owner) prepare_local();
                unsigned size = unsigned(next_history.size());
                phase = MPI_Wtime();
                mpi_check(MPI_Bcast(&size, 1, MPI_UNSIGNED, owner, comm));
                metrics.window_mpi_s += MPI_Wtime() - phase;
                if (size > HISTORY) throw std::runtime_error("invalid propagated dictionary size");
                next_history.resize(size);
                phase = MPI_Wtime();
                mpi_check(MPI_Bcast(next_history.data(), size, MPI_BYTE, owner, comm));
                metrics.window_mpi_s += MPI_Wtime() - phase;
            }
            ++metrics.serial_history_epochs;
        }
        metrics.window_local_s += MPI_Wtime() - window_begin - (metrics.window_mpi_s - window_mpi_before);
        phase = MPI_Wtime();
        if (local_bytes > config.decoded_bytes) throw std::runtime_error("decoded batch exceeds budget");
        // RESOLVE writes every byte before publication; a host zero-fill would
        // unnecessarily write the entire decoded stream a second time.
        bytes.reset(local_bytes ? new unsigned char[local_bytes] : nullptr);
        bytes_size = local_bytes;
        lines.resize(resolved.size());
        size_t offset = 0;
        unsigned char empty_output = 0;
        for (size_t i = 0; i < resolved.size(); ++i) {
            auto &job = resolved[i];
            lines[i].resize((job.output_size + PAGE - 1) / PAGE);
            job.line_counts = lines[i].data();
            job.output = job.output_size ? bytes.get() + offset : &empty_output;
            offset += job.output_size;
        }
        metrics.prepare_s += MPI_Wtime() - phase;
        launch(resolved, true);
        phase = MPI_Wtime();
        std::vector<std::array<uint64_t, 2>> digests(count), all_digests;
        if (rank == 0) all_digests.resize(size_t(ranks) * count);
        offset = 0;
        for (size_t i = 0; i < resolved.size(); ++i) {
            const auto &job = resolved[i];
            const size_t j = selected[i], id = size_t(rank) * count + j;
            if (job.status != Z_OK || job.output_size != local_meta[j][BYTES])
                throw std::runtime_error("local RESOLVE failed");
            digests[j] = {job.crc, job.output_size};
            if (job.output_size) {
                if (!extents.empty() && plan[id][1] != add(extents.front().global_offset, offset))
                    throw std::runtime_error("rank-local decoded interval is not contiguous");
                extents.push_back({job.output, job.output_size, plan[id][1], plan[id][2],
                    local_meta[j][END], lines[i].data(), lines[i].size(), job.crc});
            }
            offset += job.output_size;
            ++metrics.accepted_jobs; metrics.speculative_accepted += !job.known_start;
            metrics.decoded_bytes += job.output_size;
        }
        metrics.discarded_jobs += jobs.size() - resolved.size();
        for (size_t j = 0; j < jobs.size(); ++j)
            if (!plan[size_t(rank) * count + j][0] && local_meta[j][OK])
                metrics.discarded_output_bytes += local_meta[j][BYTES];
        metrics.commit_local_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        mpi_check(MPI_Gather(digests.data(), count * 2, MPI_UINT64_T,
            all_digests.data(), count * 2, MPI_UINT64_T, 0, comm));
        metrics.commit_mpi_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        double commit_io_before = metrics.io_s;
        // cursor, global bytes/CRC, EOF, member-ended, member CRC, true member
        // end bit, and next header fingerprint. The last two validate locally.
        std::array<uint64_t, 8> commit{};
        if (rank == 0) {
            uint32_t next_crc = crc;
            uint32_t next_member_crc = member_crc;
            for (size_t id : chain.jobs) {
                if (all_digests[id][1] != all_meta[id][BYTES] || all_digests[id][0] > UINT32_MAX)
                    throw std::runtime_error("resolve digest metadata mismatch");
                next_crc = combine.combine(next_crc, uint32_t(all_digests[id][0]), all_digests[id][1]);
                next_member_crc = combine.combine(next_member_crc, uint32_t(all_digests[id][0]), all_digests[id][1]);
            }
            commit = {chain.cursor, chain.decoded, next_crc, 0, uint64_t(chain.final), next_member_crc,
                      chain.cursor, 0};
            if (chain.final) {
                const uint64_t next = source->validate_member(chain.cursor, next_member_crc,
                                                             chain.decoded - member_begin);
                commit[3] = next == source->size();
                if (!commit[3]) {
                    uint32_t header_crc = 0;
                    commit[0] = source->header(next, header_crc);
                    commit[7] = header_crc;
                }
            }
        }
        metrics.commit_local_s += MPI_Wtime() - phase - (metrics.io_s - commit_io_before);
        // No caller sees this epoch until root validated the global commit.
        phase = MPI_Wtime();
        mpi_check(MPI_Bcast(commit.data(), commit.size(), MPI_UINT64_T, 0, comm));
        metrics.commit_mpi_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        commit_io_before = metrics.io_s;
        if (commit[4]) {
            const uint64_t next = source->validate_member(commit[6], uint32_t(commit[5]), commit[1] - member_begin);
            if (bool(commit[3]) != (next == source->size())) throw std::runtime_error("member EOF metadata mismatch");
            if (!commit[3]) {
                uint32_t header_crc = 0;
                if (source->header(next, header_crc) != commit[0] || header_crc != commit[7])
                    throw std::runtime_error("next member header metadata mismatch");
            }
            next_history.clear();
        }
        metrics.commit_local_s += MPI_Wtime() - phase - (metrics.io_s - commit_io_before);
        // All ranks finish local validation before any rank reports success.
        phase = MPI_Wtime();
        mpi_check(MPI_Barrier(comm));
        metrics.commit_barrier_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        cursor = commit[0]; decoded = commit[1]; crc = uint32_t(commit[2]); eof = commit[3] != 0;
        member_crc = commit[4] ? 0 : uint32_t(commit[5]);
        if (commit[4]) member_begin = decoded;
        history.swap(next_history);
        out.bytes = bytes.get(); out.size = bytes_size;
        out.decoded_offset = extents.empty() ? decoded : extents.front().global_offset;
        out.extents = extents.empty() ? nullptr : extents.data(); out.count = extents.size();
        out.global_decoded_end = decoded; out.stream_crc32 = crc; out.eof = eof;
        ++epoch; ++metrics.epochs;
        metrics.commit_local_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        release(symbols); release(tails);
        metrics.cleanup_s += MPI_Wtime() - phase;
        phase = MPI_Wtime();
        if (!eof) {
            const InputWindow next = input_window();
            const uint64_t retained = add(add(fixed_bytes, bytes_size), input_size);
            const uint64_t request_bytes = add(next.bytes, GUARD);
            if (next.bytes && retained < config.memory_bytes &&
                request_bytes <= config.memory_bytes - retained) {
                source->prefetch(next.base, next.bytes, input.get(), input_base,
                                 input_size ? input_size - GUARD : 0);
                metrics.peak_buffer_bytes = std::max(metrics.peak_buffer_bytes, add(retained, request_bytes));
            }
        }
        input.reset(); input_size = 0;
        metrics.prepare_s += MPI_Wtime() - phase;
        return 1;
    }
};

extern "C" void swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *c, unsigned workers) {
    if (!c) return;
    *c = swbwa_mpi_decoder_config{};
    c->chunk_bytes = 256U << 10; c->lookahead_bytes = 1U << 20;
    const size_t budget_workers = std::max(64U, std::min(workers, 384U));
    c->memory_bytes = budget_workers * (8U << 20);
    c->decoded_bytes = budget_workers * (2U << 20);
    c->initial_output_bytes = 2U << 20;
    c->jobs_per_rank = workers; c->workers = workers;
}

extern "C" int swbwa_mpi_decoder_open(MPI_Comm comm, const char *path,
    const swbwa_mpi_decoder_config *config, swbwa_mpi_decoder **out, char *error, size_t size) {
    if (out) *out = nullptr;
    if (error && size) error[0] = 0;
    try {
        if (!out) throw std::runtime_error("null decoder output");
        std::unique_ptr<swbwa_mpi_decoder> decoder(new swbwa_mpi_decoder);
        decoder->init(comm, path, config);
        *out = decoder.release();
        return 0;
    } catch (const std::exception &e) {
        // Do not free the duplicated communicator here: a peer may be blocked
        // in init. The required caller-side Abort reclaims MPI resources.
        if (error && size) std::snprintf(error, size, "%s", e.what());
        return -1;
    }
}

extern "C" int swbwa_mpi_decoder_next(swbwa_mpi_decoder *decoder, swbwa_mpi_decoded_batch *batch) {
    if (batch) *batch = swbwa_mpi_decoded_batch{};
    if (!decoder) return -1;
    try {
        if (!batch) throw std::runtime_error("null batch output");
        return decoder->next(*batch);
    } catch (const std::exception &e) {
        decoder->failed = true;
        std::snprintf(decoder->error, sizeof(decoder->error), "%s", e.what());
        if (batch) *batch = swbwa_mpi_decoded_batch{};
        return -1;
    }
}

extern "C" const char *swbwa_mpi_decoder_error(const swbwa_mpi_decoder *decoder) {
    return decoder ? decoder->error : "null MPI gzip decoder";
}

extern "C" int swbwa_mpi_decoder_get_metrics(const swbwa_mpi_decoder *decoder,
                                             swbwa_mpi_decoder_metrics *metrics) {
    if (!decoder || !metrics) return -1;
    *metrics = decoder->metrics;
    return 0;
}

extern "C" int swbwa_mpi_decoder_close(swbwa_mpi_decoder *decoder) {
    if (!decoder) return 0;
    try {
        if (decoder->failed) throw std::runtime_error("decoder failed; caller must MPI_Abort");
        if (decoder->comm != MPI_COMM_NULL) mpi_check(MPI_Comm_free(&decoder->comm));
        delete decoder;
        return 0;
    } catch (const std::exception &e) {
        decoder->failed = true;
        std::snprintf(decoder->error, sizeof(decoder->error), "%s", e.what());
        return -1;
    }
}
