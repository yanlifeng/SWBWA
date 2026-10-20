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
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

static std::thread::id owner;
extern "C" void __real_swbwa_codec_launch(swbwa_codec_batch *);
static void launch(void *p)
{
    assert(std::this_thread::get_id() == owner);
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

static void *process(void *context, int step, void *item)
{
    Fixture &f = *static_cast<Fixture *>(context);
    Batch *b = static_cast<Batch *>(item);
    if (step == 0) {
        b = new Batch;
        const int status = swbwa_mpi_gzip_next(f.reader, &b->input);
        assert(status >= 0);
        if (!status) { delete b; return NULL; }
        return b;
    }
    if (step == 1) {
        assert(std::this_thread::get_id() == owner);
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
            if (b->input.chunk_id % 5 == 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            assert(swbwa_output_write_chunk(b->input.chunk_id, b->input.first_record,
                b->input.first_record + b->input.records, b->seqs.data(), NULL,
                int(b->seqs.size()), f.mates == 2) == 0);
            assert(!swbwa_mpi_gzip_complete(f.reader, b->input.chunk_id));
        }
        delete b;
        return NULL;
    }
    assert(step == 3);
    assert(!swbwa_output_flush());
    return NULL;
}

int main(int argc, char **argv)
{
    assert(argc == 5);
    assert(!swbwa_mpi_init(&argc, &argv));
    owner = std::this_thread::get_id();
    assert(!swbwa_host_workers_init());
    swbwa_codec_init();
    Fixture f = {NULL, std::strcmp(argv[2], "-") ? 2 : 1, swbwa_mpi_rank()};
    f.reader = swbwa_mpi_gzip_open(argv[1], f.mates == 2 ? argv[2] : NULL, 2048, 64, 1);
    assert(f.reader);
    assert(!swbwa_output_open_online(argv[3], 1));
    assert(!swbwa_mpi_progress_thread_start(0));
    if (std::strcmp(argv[4], "queue") == 0) kt_pipeline_queue(3, process, &f, 3);
    else kt_pipeline_single(1, process, &f, 3);
    assert(!swbwa_output_ordered_finish(swbwa_mpi_gzip_chunk_count(f.reader)));
    assert(!swbwa_mpi_gzip_report(f.reader, stderr));
    swbwa_mpi_gzip_close(f.reader);
    assert(!swbwa_output_close());
    assert(!swbwa_mpi_progress_thread_stop());
    swbwa_host_workers_destroy(); swbwa_codec_finish();
    swbwa_mpi_finalize();
    return 0;
}
