#include <cassert>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <vector>
#include <zlib.h>
#include "codec_job.h"
#include "rapid_state.h"
#include "swbwa_host_workers.h"

#if !defined(SWBWA_CODEC_MPE) || SWBWA_CODEC_LDM != 0 || defined(SWBWA_CODEC_CPE)
#error "Test requires SWBWA_CODEC_MPE=1 and SWBWA_CODEC_LDM=0 without SWBWA_CODEC_CPE"
#endif

extern "C" int swbwa_codec_mpe_decode(swbwa_codec_batch *);
extern "C" void swbwa_codec_execute_mpe(swbwa_codec_job *);

struct alignas(256) StateProbe {
    static bool fail;
    static unsigned destroyed;
    unsigned char payload[65728];
    StateProbe() {
        assert(reinterpret_cast<uintptr_t>(this) % alignof(StateProbe) == 0);
        if (fail) throw 17;
        payload[0] = 42;
    }
    ~StateProbe() { ++destroyed; }
};
bool StateProbe::fail = false;
unsigned StateProbe::destroyed = 0;

static void test_state_lifetime()
{
    for (unsigned i = 0; i < 8; ++i) {
        {
            SwbwaRapidState<StateProbe> state;
            assert(!state.in_ldm() && state->payload[0] == 42);
        }
        assert(StateProbe::destroyed == i + 1);
        StateProbe::fail = true;
        try { SwbwaRapidState<StateProbe> state; assert(false); }
        catch (int value) { assert(value == 17); }
        assert(StateProbe::destroyed == i + 1);
        StateProbe::fail = false;
    }
    std::puts("MPE_STATE_LIFETIME_PASS alignment=256 construction_failure_cleanup=LSan");
}

int main()
{
    test_state_lifetime();
    constexpr size_t bytes = 512 * 1024, capacity = bytes + 65536;
    constexpr unsigned count = 24, iterations = 8;
    std::vector<unsigned char> raw(bytes);
    uint32_t state = 17;
    for (size_t i = 0; i < bytes; ++i) {
        state = state * 1664525U + 1013904223U;
        raw[i] = i % 151 == 150 ? '\n' : "ACGT"[(state >> 28) & 3];
    }
    std::vector<unsigned char> input(compressBound(bytes) + 131072);
    z_stream stream{};
    assert(deflateInit2(&stream, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK);
    stream.next_in = raw.data(); stream.avail_in = raw.size();
    stream.next_out = input.data(); stream.avail_out = input.size();
    assert(deflate(&stream, Z_FINISH) == Z_STREAM_END);
    size_t compressed = stream.total_out;
    assert(deflateEnd(&stream) == Z_OK);
    input.resize(compressed + 131072);
    std::vector<std::vector<uint16_t>> symbols(count, std::vector<uint16_t>(capacity));
    std::vector<std::vector<unsigned char>> output(count, std::vector<unsigned char>(capacity));
    std::vector<std::vector<uint32_t>> lines(count, std::vector<uint32_t>((capacity + 4095) / 4096));
    std::vector<swbwa_codec_job> jobs(count);
    std::vector<unsigned char> dictionary(32768, 'T');
    swbwa_codec_batch batch{};
    batch.jobs = jobs.data(); batch.count = count; batch.workers = 384;
    assert(swbwa_host_workers_init() == 0);
    for (unsigned iteration = 0; iteration < iterations; ++iteration) {
        for (unsigned i = 0; i < count; ++i) {
            auto &job = jobs[i];
            job = {};
            job.input = input.data(); job.input_size = input.size();
            job.logical_size = compressed; job.stop_bit = compressed * 8;
            job.output_capacity = capacity; job.known_start = 1;
            job.operation = SWBWA_CODEC_DISCOVER;
            job.symbols = symbols[i].data(); job.output = output[i].data();
            job.line_counts = lines[i].data();
            std::memset(job.output, 0, bytes);
        }
        assert(swbwa_codec_mpe_decode(&batch));
        for (const auto &job : jobs) {
            assert(job.status == Z_OK && job.output_size == bytes && job.final_block);
            assert(!job.state_in_ldm && job.state_bytes > 0);
            for (size_t i = 0; i < bytes; ++i) assert(job.symbols[i] == raw[i]);
        }
        for (auto &job : jobs) job.operation = SWBWA_CODEC_RESOLVE;
        assert(swbwa_codec_mpe_decode(&batch));
        for (const auto &job : jobs) {
            assert(job.status == Z_OK && !std::memcmp(job.output, raw.data(), bytes));
            assert(job.crc == crc32(0, raw.data(), raw.size()));
            size_t total = 0;
            for (size_t p = 0; p < (bytes + 4095) / 4096; ++p) total += job.line_counts[p];
            assert(total == bytes / 151);
        }
        // Resolve unresolved history references, not just literal-only input.
        for (auto &job : jobs) {
            job.dictionary = dictionary.data(); job.dictionary_size = dictionary.size();
            job.symbols[0] = 65535;
        }
        assert(swbwa_codec_mpe_decode(&batch));
        for (const auto &job : jobs) assert(job.status == Z_OK && job.output[0] == 'T');
    }
    jobs.back().operation = SWBWA_CODEC_GZIP;
    assert(!swbwa_codec_mpe_decode(&batch));
    jobs.back().operation = SWBWA_CODEC_DEFLATE;
    assert(!swbwa_codec_mpe_decode(&batch));
    swbwa_host_workers_destroy();
    std::printf("MPE_DECODE_PASS iterations=%u jobs_per_batch=%u\n", iterations, count);
}
