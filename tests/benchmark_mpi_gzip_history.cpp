// Decompression only: production decoder/CRC validation, no FASTQ or alignment.
#include "mpi_decoder.h"
#include "codec_job.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#ifdef SWBWA_HISTORY_APPLICATION_BACKEND
extern "C" {
#include "swbwa_mpi.h"
void swbwa_benchmark_runtime_init(void);
void swbwa_benchmark_runtime_finish(void);
}
#endif

extern "C" int swbwa_benchmark_cache_open(const char *) __attribute__((weak));
extern "C" void swbwa_benchmark_cache_freeze(int) __attribute__((weak));

static void checked(int status) {
    if (status != MPI_SUCCESS) throw std::runtime_error("benchmark MPI failure");
}

int main(int argc, char **argv) {
    int provided, rank, ranks;
#ifdef SWBWA_HISTORY_APPLICATION_BACKEND
    if (swbwa_mpi_init(&argc, &argv)) return 1;
    MPI_Query_thread(&provided);
#else
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try {
        if (argc < 2 || argc > 4) throw std::runtime_error("usage: benchmark FILE.gz [repeats=3] [workers=64]");
        const int repeats = argc > 2 ? std::atoi(argv[2]) : 3;
        const int workers = argc > 3 ? std::atoi(argv[3]) : 64;
        if (repeats < 1 || repeats > 10 || workers < 1 || workers > 384 || provided < MPI_THREAD_FUNNELED)
            throw std::runtime_error("invalid benchmark configuration or MPI thread level");
#ifdef SWBWA_HISTORY_APPLICATION_BACKEND
        swbwa_benchmark_runtime_init();
#endif
        swbwa_codec_init();
        const bool cached = std::getenv("SWBWA_HISTORY_BENCH_CACHE") != nullptr;
        if (cached && (!swbwa_benchmark_cache_open || !swbwa_benchmark_cache_freeze ||
                       swbwa_benchmark_cache_open(argv[1])))
            throw std::runtime_error("benchmark compressed cache unavailable");
        uint64_t expected_size = 0;
        uint32_t expected_crc = 0;
        for (int repeat = cached ? -1 : 0; repeat < repeats; ++repeat)
        for (int order = 0; order < (repeat < 0 ? 1 : 3); ++order) {
            const int mode = repeat < 0 ? 0 : (repeat + order) % 3;
            char setting[2] = {char('0' + mode), 0};
            if (setenv("SWBWA_GZIP_PARALLEL_HISTORY", setting, 1)) throw std::runtime_error("setenv failed");
            if (!rank) {
                std::printf("HISTORY_BEGIN repeat=%d mode=%d ranks=%d workers=%d file=%s\n",
                            repeat, mode, ranks, workers, argv[1]);
                std::fflush(stdout);
            }
            swbwa_mpi_decoder_config config;
            swbwa_mpi_decoder_config_init(&config, unsigned(workers));
            swbwa_mpi_decoder *decoder = nullptr;
            char error[512];
            checked(MPI_Barrier(MPI_COMM_WORLD));
            if (swbwa_mpi_decoder_open(MPI_COMM_WORLD, argv[1], &config, &decoder, error, sizeof(error)))
                throw std::runtime_error(error);
            checked(MPI_Barrier(MPI_COMM_WORLD));
            const double begin = MPI_Wtime();
            swbwa_mpi_decoded_batch batch{};
            uint64_t size = 0;
            uint32_t crc = 0;
            bool eof = false;
            for (;;) {
                const int result = swbwa_mpi_decoder_next(decoder, &batch);
                if (result < 0) throw std::runtime_error(swbwa_mpi_decoder_error(decoder));
                if (!result) break;
                size = batch.global_decoded_end; crc = batch.stream_crc32; eof = batch.eof != 0;
            }
            const double wall = MPI_Wtime() - begin;
            if (!eof) throw std::runtime_error("missing trailer-validated EOF");
            if (repeat == (cached ? -1 : 0) && !order) { expected_size = size; expected_crc = crc; }
            if (size != expected_size || crc != expected_crc) throw std::runtime_error("history modes disagree");
            swbwa_mpi_decoder_metrics metrics{};
            if (swbwa_mpi_decoder_get_metrics(decoder, &metrics)) throw std::runtime_error("metrics failed");
            if (swbwa_mpi_decoder_close(decoder)) throw std::runtime_error(swbwa_mpi_decoder_error(decoder));
            uint64_t sum = 0;
            checked(MPI_Allreduce(&metrics.decoded_bytes, &sum, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD));
            if (sum != size) throw std::runtime_error("decoded byte accounting mismatch");
            if (repeat < 0) {
                swbwa_benchmark_cache_freeze(rank);
                checked(MPI_Barrier(MPI_COMM_WORLD));
                if (!rank) std::puts("HISTORY_WARMUP_PASS timed_reads=memory_replay");
                continue;
            }
            std::vector<swbwa_mpi_decoder_metrics> all(rank ? 0 : ranks);
            std::vector<double> walls(rank ? 0 : ranks);
            checked(MPI_Gather(&metrics, sizeof(metrics), MPI_BYTE, all.data(), sizeof(metrics), MPI_BYTE, 0, MPI_COMM_WORLD));
            checked(MPI_Gather(&wall, 1, MPI_DOUBLE, walls.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
            if (!rank) {
                std::printf("HISTORY_RESULT repeat=%d mode=%d ranks=%d bytes=%llu crc32=%08x wall_max=%.9f PASS\n",
                    repeat, mode, ranks, (unsigned long long)size, unsigned(crc), *std::max_element(walls.begin(), walls.end()));
                for (int r = 0; r < ranks; ++r) {
                    const auto &m = all[r];
                    std::printf("HISTORY_RANK repeat=%d mode=%d rank=%d wall=%.9f", repeat, mode, r, walls[r]);
#define SECONDS(name) std::printf(" " #name "=%.9f", m.name)
#define COUNT(name) std::printf(" " #name "=%llu", (unsigned long long)m.name)
                    SECONDS(open_s); SECONDS(next_s); SECONDS(io_s); SECONDS(open_io_s);
                    SECONDS(source_io_s); SECONDS(prefetch_io_s); SECONDS(discover_wall_s); SECONDS(resolve_wall_s);
                    SECONDS(discover_run_s); SECONDS(resolve_run_s); SECONDS(metadata_local_s); SECONDS(metadata_mpi_s);
                    SECONDS(history_wall_s); SECONDS(window_local_s); SECONDS(window_mpi_s);
                    SECONDS(history_summary_s); SECONDS(history_exchange_s); SECONDS(history_linear_apply_s);
                    SECONDS(history_job_dict_s); SECONDS(history_tree_compose_s); SECONDS(history_tree_apply_s);
                    SECONDS(history_tree_mpi_s); SECONDS(history_final_sync_s);
                    SECONDS(commit_local_s); SECONDS(commit_mpi_s); SECONDS(commit_barrier_s); SECONDS(retry_s);
                    COUNT(epochs); COUNT(initial_batch_jobs); COUNT(max_discover_jobs); COUNT(max_resolve_jobs);
                    COUNT(discover_jobs); COUNT(accepted_jobs); COUNT(capacity_retries); COUNT(decoded_bytes);
                    COUNT(source_read_bytes); COUNT(known_tail_jobs); COUNT(repaired_tail_jobs);
                    COUNT(history_scratch_peak_bytes); COUNT(history_tree_sent_bytes); COUNT(history_tree_messages);
                    COUNT(history_tree_compose_calls); COUNT(history_tree_apply_calls); COUNT(peak_buffer_bytes);
#undef COUNT
#undef SECONDS
                    std::printf("\n");
                }
                std::fflush(stdout);
            }
        }
        swbwa_codec_finish();
#ifdef SWBWA_HISTORY_APPLICATION_BACKEND
        swbwa_benchmark_runtime_finish();
#endif
        if (!rank) std::puts("HISTORY_BENCHMARK_PASS");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "HISTORY_ERROR rank=%d %s\n", rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
#ifdef SWBWA_HISTORY_APPLICATION_BACKEND
    swbwa_mpi_finalize();
#else
    MPI_Finalize();
#endif
    return 0;
}
