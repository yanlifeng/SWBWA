#include "swbwa_mpi_gzip.h"
#include "mpi_decoder.h"
#include "codec_job.h"
#include "fastq_frame.h"
#include "swbwa_host_workers.h"
#include <mpi.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

static uint64_t decode_calls;
extern "C" void __real_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *);
extern "C" void __wrap_swbwa_mpi_decoder_config_init(swbwa_mpi_decoder_config *c)
{
    __real_swbwa_mpi_decoder_config_init(c);
    if (!std::getenv("FASTQ_TEST_DEFAULT_CONFIG")) {
        c->chunk_bytes = 4096; c->jobs_per_rank = 1; c->initial_output_bytes = 65536;
        c->decoded_bytes = 8U << 20; c->memory_bytes = 64U << 20;
    }
    if (const char *n = std::getenv("FASTQ_TEST_LOOKAHEAD"))
        c->lookahead_bytes = size_t(std::strtoull(n, NULL, 10));
}
extern "C" int __real_swbwa_mpi_decoder_next(swbwa_mpi_decoder *, swbwa_mpi_decoded_batch *);
extern "C" int __wrap_swbwa_mpi_decoder_next(swbwa_mpi_decoder *d, swbwa_mpi_decoded_batch *b)
{
    ++decode_calls; return __real_swbwa_mpi_decoder_next(d, b);
}
static std::string read(const char *path)
{
    std::ifstream file(path, std::ios::binary); assert(file);
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}
static void validate_partitions(const swbwa_mpi_gzip_batch &b, bool paired)
{
    assert(b.partitions && b.partition_count);
    size_t previous[2] = {0, 0};
    int64_t ordinal = 0;
    for (unsigned c = 0; c < b.partition_count; ++c) {
        const swbwa_fastq_partition_t &part = b.partitions[c];
        assert(part.first_record == ordinal && part.records >= 0);
        for (int m = 0; m < 1 + int(paired); ++m)
            assert(part.begin[m] == int64_t(previous[m]) && part.end[m] >= part.begin[m] && part.end[m] <= b.size[m]);
        size_t positions[2] = {previous[0], previous[1]};
        for (int64_t i = 0; i < part.records; ++i) {
            swbwa_fastq_record_view record[2];
            for (int m = 0; m < 1 + int(paired); ++m) {
                if (swbwa_fastq_record_next(b.data[m], size_t(part.end[m]), &positions[m], &record[m]) != 1) {
                    std::fprintf(stderr, "invalid FASTQ in CPE validation\n");
                    MPI_Abort(MPI_COMM_WORLD, 1); std::abort();
                }
            }
            if (paired && (record[0].name_size != record[1].name_size ||
                std::memcmp(record[0].name, record[1].name, record[0].name_size))) {
                std::fprintf(stderr, "paired FASTQ name token mismatch in CPE validation\n");
                MPI_Abort(MPI_COMM_WORLD, 1); std::abort();
            }
        }
        for (int m = 0; m < 1 + int(paired); ++m) {
            assert(positions[m] == size_t(part.end[m])); previous[m] = positions[m];
        }
        ordinal += part.records;
    }
    assert(ordinal == b.records && previous[0] == size_t(b.size[0]));
    if (paired) assert(previous[1] == size_t(b.size[1]));
}
static void run(int argc, char **argv)
{
    assert(argc == 9);
    int rank, ranks; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const bool paired = std::strcmp(argv[2], "-") != 0;
    const std::string mode = argv[8];
    const bool rejection = mode == "reject" || mode == "early-close" || mode == "config-mismatch";
    std::string reference[2];
    std::vector<swbwa_fastq::Record> records[2];
    if (!rejection) {
        for (int m = 0; m < 1 + int(paired); ++m) {
            reference[m] = read(argv[3 + m]);
            swbwa_fastq::LineOffsets lines;
            lines.scan(reinterpret_cast<const unsigned char *>(reference[m].data()), reference[m].size());
            records[m] = swbwa_fastq::parse(reinterpret_cast<const unsigned char *>(reference[m].data()), lines, 0, reference[m].size());
        }
    }
    size_t target = size_t(std::strtoull(argv[5], NULL, 10));
    if (mode == "config-mismatch" && rank) ++target;
    swbwa_mpi_gzip_reader *reader = swbwa_mpi_gzip_open(argv[1], paired ? argv[2] : NULL,
        target, unsigned(std::atoi(argv[6])), std::atoi(argv[7]));
    assert(reader && decode_calls == 0); // open must not sample/decode.
    if (mode == "early-close") swbwa_mpi_gzip_close(reader);
    std::vector<int> coverage(records[0].size(), 0), total(records[0].size(), 0);
    uint64_t rounds = 0, local_records = 0, empty = 0;
    char *held[2] = {NULL, NULL}; size_t held_size[2] = {0, 0};
    for (;;) {
        swbwa_mpi_gzip_batch b = {};
        const int status = swbwa_mpi_gzip_next(reader, &b);
        int low = 0, high = 0;
        MPI_Allreduce(&status, &low, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&status, &high, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        assert(low == high && status >= 0);
        for (int m = 0; m < 2; ++m) {
            if (held[m]) {
                for (size_t i = 0; i < held_size[m]; ++i) assert(static_cast<unsigned char>(held[m][i]) == 0xa5);
                std::free(held[m]); held[m] = NULL;
            }
        }
        if (!status) {
            assert(!b.data[0] && !b.data[1] && !b.records);
            assert(!swbwa_mpi_gzip_next(reader, &b)); break;
        }
        assert(b.chunk_id == int64_t(rounds * ranks + rank)); ++rounds;
        if (!b.records) {
            assert(!b.data[0] && !b.data[1] && !b.size[0] && !b.size[1] && !b.partitions); ++empty; continue;
        }
        validate_partitions(b, paired);
        std::free(b.partitions);
        if (rejection) {
            std::free(b.data[0]); std::free(b.data[1]); continue;
        }
        assert(b.first_record >= 0 && b.records > 0 && uint64_t(b.first_record + b.records) <= coverage.size());
        for (int64_t i = b.first_record; i < b.first_record + b.records; ++i) ++coverage[size_t(i)];
        local_records += uint64_t(b.records);
        for (int m = 0; m < 1 + int(paired); ++m) {
            const size_t begin = records[m][size_t(b.first_record)].begin;
            const size_t end = records[m][size_t(b.first_record + b.records - 1)].end;
            assert(b.data[m] && b.size[m] == int64_t(end - begin) && b.data[m][b.size[m]] == 0);
            assert(!std::memcmp(b.data[m], reference[m].data() + begin, end - begin));
            held[m] = b.data[m]; held_size[m] = size_t(b.size[m]);
            std::memset(held[m], 0xa5, held_size[m]);
        }
        if (!paired) assert(!b.data[1] && !b.size[1]);
    }
    assert(!rejection);
    assert(!*swbwa_mpi_gzip_error(reader));
    assert(!swbwa_mpi_gzip_report(reader, stderr));
    swbwa_mpi_gzip_close(reader);
    if (!coverage.empty()) {
        MPI_Allreduce(coverage.data(), total.data(), int(total.size()), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for (int n : total) assert(n == 1);
    }
    unsigned long long own[3] = {local_records, empty, decode_calls}, sum[3] = {};
    MPI_Reduce(own, sum, 3, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank) {
        assert(sum[0] == records[0].size());
        std::printf("MPI_FASTQ_RESULT {\"records\":%llu,\"rounds\":%llu,\"empty_ready\":%llu,\"decode_calls\":%llu}\n",
                    sum[0], static_cast<unsigned long long>(rounds), sum[1], sum[2]);
    }
}
int main(int argc, char **argv)
{
    int provided = 0;
    assert(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided) == MPI_SUCCESS);
    assert(provided >= MPI_THREAD_MULTIPLE);
    swbwa_codec_init();
    assert(swbwa_host_workers_init() == 0);
    if (argc == 9 && std::strcmp(argv[8], "thread") == 0) {
        std::thread worker([&]() { run(argc, argv); }); worker.join();
    } else run(argc, argv);
    swbwa_host_workers_destroy();
    swbwa_codec_finish(); MPI_Finalize();
}
