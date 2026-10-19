#include "mpi_decoder.h"
#include "mpi_decoder_protocol.h"
#include "codec_job.h"
#include <zlib.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

static unsigned largest_launch = 0;
static bool fail_comm_free = false;
static int stat_time_skew = 0;
static unsigned stat_time_skew_calls = 0;
static int skew_stat_time(int result, struct stat64 *info) {
    if (!result && stat_time_skew) {
        info->st_mtim.tv_sec -= stat_time_skew;
        info->st_ctim.tv_sec -= stat_time_skew;
        ++stat_time_skew_calls;
    }
    return result;
}
extern "C" int __real_fstat64(int, struct stat64 *);
extern "C" int __wrap_fstat64(int fd, struct stat64 *info) {
    return skew_stat_time(__real_fstat64(fd, info), info);
}
#if defined(__GLIBC__)
#if !__GLIBC_PREREQ(2, 33)
// Older glibc headers inline fstat64 through the versioned entry point.
extern "C" int __real___fxstat64(int, int, struct stat64 *);
extern "C" int __wrap___fxstat64(int version, int fd, struct stat64 *info) {
    return skew_stat_time(__real___fxstat64(version, fd, info), info);
}
#endif
#endif
extern "C" int __real_MPI_Comm_free(MPI_Comm *);
extern "C" int __wrap_MPI_Comm_free(MPI_Comm *comm) {
    if (fail_comm_free) { fail_comm_free = false; return MPI_ERR_OTHER; }
    return __real_MPI_Comm_free(comm);
}
extern "C" void __real_swbwa_codec_launch(swbwa_codec_batch *);
extern "C" void __wrap_swbwa_codec_launch(swbwa_codec_batch *batch) {
    assert(batch->workers >= 1 && batch->workers <= 384);
    largest_launch = std::max(largest_launch, batch->count);
    __real_swbwa_codec_launch(batch);
}

static void require(bool ok, const char *what) {
    if (!ok) throw std::runtime_error(what);
}

static void unit() {
    using namespace swbwa_mpi_gzip;
    std::vector<Meta> metas(6);
    for (unsigned i = 0; i < 6; ++i) {
        const unsigned id = (4 + i) % 6;
        auto &m = metas[id];
        m[PRESENT] = m[OK] = 1; m[KNOWN] = i == 0;
        m[BASE] = m[START] = m[START_MAX] = 80 + i * 10;
        m[END] = m[LIMIT] = 90 + i * 10;
        m[BYTES] = i == 1 ? 0 : 7;
    }
    metas[3][FINAL] = 1;
    std::vector<Plan> plan;
    auto chain = select(metas, plan, 3, 2, 2, 80, 100, 200, 65536);
    assert(chain.jobs == std::vector<size_t>({4, 5, 0, 1, 2, 3}));
    assert(chain.final && chain.cursor == 140 && chain.decoded == 135);
    assert(plan[0][1] == 107 && plan[4][1] == 100);
    metas[5][OK] = 0;
    metas[4][END] = metas[4][LIMIT] = 100;
    metas[0][START] = 98; metas[0][BASE] = 98; metas[0][START_MAX] = 102;
    chain = select(metas, plan, 3, 2, 2, 80, 0, 200, 65536);
    assert(chain.jobs.size() == 5 && plan[0][2] == 100);
    Window previous{10, 11, 12};
    assert(advance(previous, nullptr, 0) == previous);
    uint16_t tail[] = {65533, 65535, 65};
    assert(advance(previous, tail, 3) == Window({10, 11, 12, 10, 12, 65}));
    bool failed = false;
    try { uint16_t bad = 256; advance(previous, &bad, 1); }
    catch (const std::exception &) { failed = true; }
    assert(failed);
    failed = false;
    try { add(UINT64_MAX, 1); } catch (const std::exception &) { failed = true; }
    assert(failed);
    std::puts("MPI_GZIP_PROTOCOL_PASS");
}

static void verify(const swbwa_mpi_decoded_batch &batch, FILE *reference) {
    require(batch.size == 0 || (batch.bytes && batch.count), "missing contiguous buffer");
    size_t at = 0;
    for (size_t i = 0; i < batch.count; ++i) {
        const auto &e = batch.extents[i];
        require(e.data == batch.bytes + at && e.global_offset == batch.decoded_offset + at,
                "non-contiguous local extent");
        require(e.size && e.line_count_pages == (e.size + 4095) / 4096, "bad newline index size");
        for (size_t page = 0; page < e.line_count_pages; ++page) {
            const size_t begin = page * 4096, end = std::min(e.size, begin + 4096);
            uint32_t count = 0;
            for (size_t j = begin; j < end; ++j) count += e.data[j] == '\n';
            require(count == e.line_counts[page], "newline page mismatch");
        }
        require(crc32(0, e.data, e.size) == e.crc32, "local CRC mismatch");
        at += e.size;
    }
    require(at == batch.size, "incomplete local extent coverage");
    if (!batch.size) return;
    require(fseeko(reference, off_t(batch.decoded_offset), SEEK_SET) == 0, "seek reference");
    unsigned char expected[8192];
    for (at = 0; at < batch.size;) {
        size_t n = std::min(sizeof(expected), batch.size - at);
        require(fread(expected, 1, n, reference) == n, "short reference");
        require(std::memcmp(expected, batch.bytes + at, n) == 0, "decoded byte mismatch");
        at += n;
    }
}

int main(int argc, char **argv) {
    if (argc == 2 && !std::strcmp(argv[1], "--unit")) { unit(); return 0; }
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int world_rank = 0, n = 0, rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); MPI_Comm_size(MPI_COMM_WORLD, &n);
    try {
        require(argc == 11, "source reference jobs capacity chunk lookahead memory decoded mode");
        swbwa_mpi_decoder_config cfg;
        swbwa_mpi_decoder_config_init(&cfg);
        cfg.jobs_per_rank = unsigned(std::strtoul(argv[3], nullptr, 10));
        cfg.initial_output_bytes = std::strtoull(argv[4], nullptr, 10);
        cfg.chunk_bytes = std::strtoull(argv[5], nullptr, 10);
        cfg.lookahead_bytes = std::strtoull(argv[6], nullptr, 10);
        cfg.memory_bytes = std::strtoull(argv[7], nullptr, 10);
        cfg.decoded_bytes = std::strtoull(argv[8], nullptr, 10);
        const std::string mode = argv[9];
        struct stat original{};
        if (mode == "touch") require(stat(argv[1], &original) == 0, "initial timestamp stat");
        if (mode == "stat-time-skew") stat_time_skew = 60 * (world_rank + 1);
        cfg.workers = unsigned(std::strtoul(argv[10], nullptr, 10));
        MPI_Comm parent = MPI_COMM_NULL;
        MPI_Comm_split(MPI_COMM_WORLD, 0, mode == "reverse" ? n - world_rank : world_rank, &parent);
        MPI_Comm_rank(parent, &rank);
        if (mode == "config-mismatch" && rank == n - 1) ++cfg.chunk_bytes;
        if (mode == "invalid-local" && rank == n - 1) cfg.workers = 0;
        std::string path = argv[1];
        if (mode == "path-mismatch" && rank == n - 1) path += ".missing";
        swbwa_mpi_decoder *decoder = nullptr, *second = nullptr;
        char error[512];
        if (swbwa_mpi_decoder_open(parent, path.c_str(), &cfg, &decoder, error, sizeof(error)))
            throw std::runtime_error(error);
        if (mode == "dual" && swbwa_mpi_decoder_open(parent, path.c_str(), &cfg, &second, error, sizeof(error)))
            throw std::runtime_error(error);
        FILE *reference = std::fopen(argv[2], "rb");
        require(reference != nullptr, "open reference");
        uint64_t end = 0, epoch = 0;
        uint32_t crc = 0;
        bool finished = false;
        bool timestamp_changed = false;
        for (;;) {
            // Keep an application collective pending on WORLD while decoding
            // exclusively on the private duplicated input communicator.
            int value = 1, sum = 0;
            MPI_Request pending;
            MPI_Iallreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &pending);
            swbwa_mpi_decoded_batch batch{};
            int status = swbwa_mpi_decoder_next(decoder, &batch);
            if (status < 0) throw std::runtime_error(swbwa_mpi_decoder_error(decoder));
            if (second) {
                swbwa_mpi_decoded_batch other{};
                int other_status = swbwa_mpi_decoder_next(second, &other);
                if (other_status < 0) throw std::runtime_error(swbwa_mpi_decoder_error(second));
                require(other_status == status && batch.size == other.size &&
                    batch.decoded_offset == other.decoded_offset && batch.eof == other.eof &&
                    (!batch.size || !std::memcmp(batch.bytes, other.bytes, batch.size)), "dual reader mismatch");
                verify(other, reference);
            }
            MPI_Wait(&pending, MPI_STATUS_IGNORE);
            require(sum == n, "application communicator collision");
            uint64_t state[5] = {uint64_t(status), uint64_t(batch.eof), batch.epoch,
                                 batch.global_decoded_end, batch.stream_crc32};
            uint64_t low[5], high[5];
            MPI_Allreduce(state, low, 5, MPI_UINT64_T, MPI_MIN, parent);
            MPI_Allreduce(state, high, 5, MPI_UINT64_T, MPI_MAX, parent);
            require(!std::memcmp(low, high, sizeof(low)), "ranks disagree on committed state");
            if (!status) {
                require(finished && batch.eof && !batch.bytes && !batch.count, "bad repeated EOF");
                swbwa_mpi_decoded_batch again{};
                require(swbwa_mpi_decoder_next(decoder, &again) == 0 && again.eof, "unstable EOF");
                break;
            }
            require(!finished && batch.epoch == epoch && batch.global_begin == end, "bad epoch ordering");
            verify(batch, reference);
            uint64_t interval[2] = {batch.decoded_offset, batch.size};
            std::vector<uint64_t> intervals(size_t(n) * 2);
            MPI_Allgather(interval, 2, MPI_UINT64_T, intervals.data(), 2, MPI_UINT64_T, parent);
            uint64_t cursor = end;
            for (int i = 0; i < n; ++i) {
                const int owner = int((epoch + i) % n);
                if (!intervals[owner * 2 + 1]) continue;
                require(intervals[owner * 2] == cursor, "missing/duplicated/interleaved bytes");
                cursor += intervals[owner * 2 + 1];
            }
            require(cursor == batch.global_decoded_end, "incomplete global byte coverage");
            end = cursor; crc = batch.stream_crc32; finished = batch.eof; ++epoch;
            if (mode == "touch" && epoch == 1 && !finished) {
                if (!rank) {
                    struct timespec times[] = {original.st_atim, original.st_mtim};
                    times[1].tv_sec += 60;
                    require(utimensat(AT_FDCWD, argv[1], times, 0) == 0, "change timestamp");
                }
                MPI_Barrier(parent);
                timestamp_changed = true;
            }
        }
        if (mode == "touch") require(timestamp_changed, "timestamp change was not exercised");
        if (mode == "stat-time-skew") require(stat_time_skew_calls > 0, "stat skew was not exercised");
        require(fseeko(reference, 0, SEEK_END) == 0 && uint64_t(ftello(reference)) == end, "wrong final length");
        std::fclose(reference);
        swbwa_mpi_decoder_metrics m{};
        require(!swbwa_mpi_decoder_get_metrics(decoder, &m), "metrics failed");
        require(m.peak_buffer_bytes <= cfg.memory_bytes && m.epochs == epoch, "metrics budget mismatch");
        require(m.discover_jobs == m.accepted_jobs + m.discarded_jobs + m.retry_jobs,
                "discover job accounting mismatch");
        require(m.discover_output_full_jobs <= m.discover_failed_jobs &&
                m.discover_failed_jobs <= m.discover_jobs, "failed job accounting mismatch");
        const double phases = m.io_s - m.open_io_s + m.prepare_s + m.cleanup_s +
            m.discover_wall_s + m.resolve_wall_s + m.metadata_local_s + m.metadata_mpi_s +
            m.window_local_s + m.window_mpi_s + m.commit_local_s + m.commit_mpi_s +
            m.commit_barrier_s + m.eof_barrier_s;
        const double tolerance = 1e-6;
        require(phases <= m.next_s + tolerance, "overlapping next phase timers");
        require(m.open_io_s <= m.open_s + tolerance && m.open_io_s <= m.io_s + tolerance &&
                m.retry_s <= m.next_s + tolerance, "inclusive timer accounting mismatch");
        require(m.discover_run_s + m.discover_wait_s <= m.discover_wall_s + tolerance &&
                m.resolve_run_s + m.resolve_wait_s <= m.resolve_wall_s + tolerance,
                "codec run/wait not contained in wall time");
        const char *timer_names[] = {"open_s", "open_io_s", "next_s", "io_s", "prepare_s", "cleanup_s",
            "discover_wall_s", "discover_run_s", "discover_wait_s", "resolve_wall_s", "resolve_run_s",
            "resolve_wait_s", "metadata_local_s", "metadata_mpi_s", "window_local_s", "window_mpi_s",
            "commit_local_s", "commit_mpi_s", "commit_barrier_s", "eof_barrier_s", "retry_s"};
        const double local_times[] = {m.open_s, m.open_io_s, m.next_s, m.io_s, m.prepare_s, m.cleanup_s,
            m.discover_wall_s, m.discover_run_s, m.discover_wait_s, m.resolve_wall_s, m.resolve_run_s,
            m.resolve_wait_s, m.metadata_local_s, m.metadata_mpi_s, m.window_local_s, m.window_mpi_s,
            m.commit_local_s, m.commit_mpi_s, m.commit_barrier_s, m.eof_barrier_s, m.retry_s};
        constexpr int timer_count = sizeof(local_times) / sizeof(local_times[0]);
        static_assert(timer_count == sizeof(timer_names) / sizeof(timer_names[0]), "timer names mismatch");
        for (double seconds : local_times)
            require(std::isfinite(seconds) && seconds >= 0, "invalid timer value");
        const char *work_names[] = {"discover_jobs", "retry_jobs", "discover_failed_jobs",
                                   "discover_output_full_jobs", "discarded_output_bytes"};
        const uint64_t local_work[] = {m.discover_jobs, m.retry_jobs, m.discover_failed_jobs,
                                      m.discover_output_full_jobs, m.discarded_output_bytes};
        constexpr int work_count = sizeof(local_work) / sizeof(local_work[0]);
        std::vector<double> times(size_t(n) * timer_count);
        std::vector<uint64_t> work(size_t(n) * work_count);
        MPI_Gather(local_times, timer_count, MPI_DOUBLE, times.data(), timer_count, MPI_DOUBLE, 0, parent);
        MPI_Gather(local_work, work_count, MPI_UINT64_T, work.data(), work_count, MPI_UINT64_T, 0, parent);
        uint64_t local[7] = {m.accepted_jobs, m.speculative_accepted, m.discarded_jobs,
            m.capacity_retries, m.decoded_bytes, largest_launch, m.peak_buffer_bytes};
        std::vector<uint64_t> stats(size_t(n) * 7);
        MPI_Gather(local, 7, MPI_UINT64_T, stats.data(), 7, MPI_UINT64_T, 0, parent);
        if (swbwa_mpi_decoder_close(second))
            throw std::runtime_error(std::string("decoder close failed: ") + swbwa_mpi_decoder_error(second));
        fail_comm_free = mode == "close-error" && rank == 0;
        if (swbwa_mpi_decoder_close(decoder))
            throw std::runtime_error(std::string("decoder close failed: ") + swbwa_mpi_decoder_error(decoder));
        if (!rank) {
            std::printf("MPI_GZIP_RESULT {\"ranks\":%d,\"epochs\":%llu,\"decoded_bytes\":%llu,\"crc32\":%u,\"stats\":[",
                n, (unsigned long long)epoch, (unsigned long long)end, crc);
            for (int r = 0; r < n; ++r) {
                std::printf("%s[", r ? "," : "");
                for (int i = 0; i < 7; ++i) std::printf("%s%llu", i ? "," : "", (unsigned long long)stats[r * 7 + i]);
                std::printf("]");
            }
            std::printf("],\"timers\":[");
            for (int r = 0; r < n; ++r) {
                std::printf("%s{", r ? "," : "");
                for (int i = 0; i < timer_count; ++i)
                    std::printf("%s\"%s\":%.9f", i ? "," : "", timer_names[i], times[r * timer_count + i]);
                std::printf("}");
            }
            std::printf("],\"work\":[");
            for (int r = 0; r < n; ++r) {
                std::printf("%s{", r ? "," : "");
                for (int i = 0; i < work_count; ++i)
                    std::printf("%s\"%s\":%llu", i ? "," : "", work_names[i],
                                (unsigned long long)work[r * work_count + i]);
                std::printf("}");
            }
            std::puts("]}");
        }
        MPI_Comm_free(&parent); MPI_Finalize(); return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "MPI_GZIP_ERROR rank=%d %s\n", world_rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 2); return 2;
    }
}
