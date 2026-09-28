#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <thread>
#include <vector>
#include "swbwa_config.h"

extern "C" void kt_pipeline_queue(int, void* (*)(void*, int, void*), void*, int);
extern "C" void kt_pipeline_wait_report(void);

struct Batch { size_t id; uint64_t* slot; };
struct Fixture {
    size_t total, produced = 0, processed = 0, written = 0;
    unsigned flushed = 0, slow;
    uint64_t slots[SWBWA_PIPELINE_BUFFER_COUNT] = {};
    std::atomic<bool> occupied[SWBWA_PIPELINE_BUFFER_COUNT];
    explicit Fixture(size_t n, unsigned mode) : total(n), slow(mode) {
        for (auto& flag : occupied) flag = false;
    }
};

static void* process(void* opaque, int step, void* data)
{
    Fixture& f = *static_cast<Fixture*>(opaque);
    Batch* batch = static_cast<Batch*>(data);
    if (step == 0) {
        if (f.produced == f.total) return nullptr;
        if (f.slow == 1 && f.produced % 31 == 0)
            std::this_thread::sleep_for(std::chrono::microseconds(50));
        return new Batch{f.produced++, nullptr};
    }
    if (step == 1) {
        assert(batch->id == f.processed++);
        size_t slot = batch->id % SWBWA_PIPELINE_BUFFER_COUNT;
        assert(!f.occupied[slot].exchange(true));
        batch->slot = &f.slots[slot];
        *batch->slot = batch->id ^ UINT64_C(0xdeadbeef);
        return batch;
    }
    if (step == 2) {
        assert(batch->id == f.written++);
        if (f.slow == 2 && batch->id % 17 == 0)
            std::this_thread::sleep_for(std::chrono::microseconds(50));
        assert(*batch->slot == (batch->id ^ UINT64_C(0xdeadbeef)));
        f.occupied[batch->id % SWBWA_PIPELINE_BUFFER_COUNT] = false;
        delete batch;
        return nullptr;
    }
    assert(step == 3 && data == nullptr);
    ++f.flushed;
    return nullptr;
}

int main()
{
    for (unsigned mode = 0; mode < 3; ++mode) {
        for (size_t count : {size_t(0), size_t(1), size_t(7), size_t(5000)}) {
            Fixture f(count, mode);
            kt_pipeline_queue(3, process, &f, 3);
            assert(f.produced == count && f.processed == count && f.written == count && f.flushed == 1);
            for (auto& flag : f.occupied) assert(!flag);
        }
    }
    kt_pipeline_wait_report();
}
