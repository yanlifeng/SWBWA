#include <iostream>
#include <vector>
#include <thread>
#include <functional>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <pthread.h>
#include <cstring>
#include <chrono>
#include <cstdio>
#include "swbwa_pipeline.h"

#include <cassert>


/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
    struct kt_for_t *t;
    long i;
} ktf_worker_t;

typedef struct kt_for_t {
    int n_threads;
    long n;
    ktf_worker_t *w;
    void (*func)(void*,long,int);
    void *data;
} kt_for_t;


void kt_for_single(int n_threads, void (*func)(void*, long, int), void* data, long n) {
    for (int i = 0; i < n; i++) {
        func(data, i, 0);
    }
}

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
    struct ktp_t *pl;
    int64_t index;
    int step;
    void *data;
} ktp_worker_t;

typedef struct ktp_t {
    void *shared;
    void *(*func)(void*, int, void*);
    int64_t index;
    int n_workers, n_steps;
    ktp_worker_t *workers;
    pthread_mutex_t mutex;
    pthread_cond_t cv;
} ktp_t;


static void ktp_worker_single(ktp_t* p) {
    int step = 0;
    void* data = 0;
    while (step < p->n_steps) {
        data = p->func(p->shared, step, step ? data : 0);
        step = (step == p->n_steps - 1 || data) ? (step + 1) % p->n_steps : p->n_steps;
    }
    p->func(p->shared, 3, 0);
}

extern "C" void kt_pipeline_single(int n_threads, void* (*func)(void*, int, void*), void* shared_data, int n_steps) {
    ktp_t aux;
    aux.n_steps = n_steps;
    aux.func = func;
    aux.shared = shared_data;
    ktp_worker_single(&aux);
}

//void ktp_worker(ktp_worker_t* w) {
//    ktp_t* p = w->pl;
//    while (w->step < p->n_steps) {
//        std::unique_lock<std::mutex> lock(p->mtx);
//        while (true) {
//            bool can_start = true;
//            for (int i = 0; i < p->n_workers; ++i) {
//                if (w == &p->workers[i]) continue;
//                if (p->workers[i].step <= w->step && p->workers[i].index < w->index) {
//                    can_start = false;
//                    break;
//                }
//            }
//            if (can_start) break;
//            p->cv.wait(lock);
//        }
//        lock.unlock();
//
//        w->data = p->func(p->shared, w->step, w->step ? w->data : nullptr);
//
//        lock.lock();
//        w->step = (w->step == p->n_steps - 1 || w->data) ? (w->step + 1) % p->n_steps : p->n_steps;
//        if (w->step == 0) w->index = p->index++;
//        p->cv.notify_all();
//        lock.unlock();
//    }
//}
//
//extern "C" void kt_pipeline(int n_threads, void* (*func)(void*, int, void*), void* shared_data, int n_steps) {
//    ktp_t aux;
//    if (n_threads < 1) n_threads = 1;
//
//    aux.n_workers = n_threads;
//    aux.n_steps = n_steps;
//    aux.func = func;
//    aux.shared = shared_data;
//    aux.index = 0;
//    aux.workers.resize(n_threads);
//
//    for (int i = 0; i < n_threads; ++i) {
//        ktp_worker_t& w = aux.workers[i];
//        w.step = 0;
//        w.pl = &aux;
//        w.data = nullptr;
//        w.index = aux.index++;
//    }
//
//    std::vector<std::thread> threads;
//    for (int i = 0; i < n_threads; ++i) {
//        threads.emplace_back(ktp_worker, &aux.workers[i]);
//    }
//
//    for (auto& t : threads) {
//        t.join();
//    }
//}



using DataType = void*;

const int read_queue_item_limit = SWBWA_PIPELINE_INPUT_QUEUE_CAPACITY;
const int write_queue_item_limit = SWBWA_PIPELINE_QUEUE_CAPACITY;
static_assert(SWBWA_PIPELINE_BUFFER_COUNT >=
              SWBWA_PIPELINE_QUEUE_CAPACITY + 1,
              "SAM buffer ring must cover the output queue and active batch");



struct QueueWait {
    uint64_t events = 0;
    double seconds = 0;
};

static struct {
    QueueWait reader_full, processor_empty, processor_full, writer_empty;
    bool enabled = false;
} pipeline_waits;

static void pipeline_thread_check(int error)
{
    if (error == 0) return;
    fprintf(stderr, "pipeline queue synchronization failed: %s\n", strerror(error));
    abort();
}

/* Both queues share a condition variable so a reader OR writer codec request
 * can wake the owner, including during output backpressure and final draining. */
class PipelineSync {
    struct Call {
        void (*function)(void *);
        void *argument;
        Call *next;
        bool done;
    };
    Call *first = nullptr, *last = nullptr;
    bool writer_done = false;

public:
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t changed = PTHREAD_COND_INITIALIZER;

    static void checked(int error)
    {
        pipeline_thread_check(error);
    }

    ~PipelineSync()
    {
        assert(!first);
        checked(pthread_cond_destroy(&changed));
        checked(pthread_mutex_destroy(&mutex));
    }
    void service()
    {
        while (first) {
            Call *request = first;
            first = request->next;
            if (!first) last = nullptr;
            checked(pthread_mutex_unlock(&mutex));
            request->function(request->argument);
            checked(pthread_mutex_lock(&mutex));
            request->done = true;
            checked(pthread_cond_broadcast(&changed));
        }
    }
    void call(void (*function)(void *), void *argument)
    {
        Call request{function, argument, nullptr, false};
        checked(pthread_mutex_lock(&mutex));
        if (last) last->next = &request;
        else first = &request;
        last = &request;
        checked(pthread_cond_broadcast(&changed));
        while (!request.done) checked(pthread_cond_wait(&changed, &mutex));
        checked(pthread_mutex_unlock(&mutex));
    }
    template<class Predicate>
    void wait(Predicate ready, QueueWait& stats, bool owner)
    {
        bool counted = false;
        for (;;) {
            if (owner) service();
            if (ready()) return;
            const auto start = std::chrono::steady_clock::now();
            checked(pthread_cond_wait(&changed, &mutex));
            if (!counted) { ++stats.events; counted = true; }
            stats.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
    }
    void finish_writer()
    {
        checked(pthread_mutex_lock(&mutex));
        writer_done = true;
        checked(pthread_cond_broadcast(&changed));
        checked(pthread_mutex_unlock(&mutex));
    }
    void drain()
    {
        QueueWait ignored;
        checked(pthread_mutex_lock(&mutex));
        wait([this] { return writer_done; }, ignored, true);
        checked(pthread_mutex_unlock(&mutex));
    }
};

/* Single producer/consumer; closure and the last item share one mutex. */
class BatchQueue {
    std::vector<DataType> items;
    size_t head = 0, tail = 0, count = 0;
    bool closed = false;
    PipelineSync& sync;

public:
    BatchQueue(size_t capacity, PipelineSync& shared) : items(capacity), sync(shared) { assert(capacity > 0); }
    BatchQueue(const BatchQueue&) = delete;
    BatchQueue& operator=(const BatchQueue&) = delete;

    void wait_for_space(QueueWait& stats, bool owner = false)
    {
        PipelineSync::checked(pthread_mutex_lock(&sync.mutex));
        sync.wait([this] { return count < items.size(); }, stats, owner);
        PipelineSync::checked(pthread_mutex_unlock(&sync.mutex));
    }

    void push(DataType item)
    {
        PipelineSync::checked(pthread_mutex_lock(&sync.mutex));
        assert(!closed && count < items.size());
        items[tail] = item;
        tail = (tail + 1) % items.size();
        ++count;
        PipelineSync::checked(pthread_cond_broadcast(&sync.changed));
        PipelineSync::checked(pthread_mutex_unlock(&sync.mutex));
    }

    bool pop(DataType& item, QueueWait& stats, bool owner = false)
    {
        PipelineSync::checked(pthread_mutex_lock(&sync.mutex));
        sync.wait([this] { return count != 0 || closed; }, stats, owner);
        if (!count) {
            PipelineSync::checked(pthread_mutex_unlock(&sync.mutex));
            return false;
        }
        item = items[head];
        head = (head + 1) % items.size();
        --count;
        PipelineSync::checked(pthread_cond_broadcast(&sync.changed));
        PipelineSync::checked(pthread_mutex_unlock(&sync.mutex));
        return true;
    }

    void close()
    {
        PipelineSync::checked(pthread_mutex_lock(&sync.mutex));
        closed = true;
        PipelineSync::checked(pthread_cond_broadcast(&sync.changed));
        PipelineSync::checked(pthread_mutex_unlock(&sync.mutex));
    }
};

/* Compiler TLS was observed corrupted after native alignment. Keep pipeline
 * handoff pointers in pthread-managed storage for both reader and writer. */
static pthread_key_t reader_processor_key;
static pthread_once_t reader_processor_key_once = PTHREAD_ONCE_INIT;

static void create_reader_processor_key()
{
    pipeline_thread_check(pthread_key_create(&reader_processor_key, nullptr));
}

static void prepare_reader_processor_key()
{
    pipeline_thread_check(pthread_once(&reader_processor_key_once, create_reader_processor_key));
}

extern "C" void swbwa_pipeline_call_on_processor(void (*call)(void *), void *argument)
{
    assert(call);
    prepare_reader_processor_key();
    PipelineSync *sync = static_cast<PipelineSync *>(pthread_getspecific(reader_processor_key));
    if (sync)
        sync->call(call, argument);
    else
        call(argument);
}

extern "C" void *swbwa_pipeline_context(void)
{
    prepare_reader_processor_key();
    return pthread_getspecific(reader_processor_key);
}

extern "C" void swbwa_pipeline_set_context(void *context)
{
    prepare_reader_processor_key();
    pipeline_thread_check(pthread_setspecific(reader_processor_key, context));
}

static void reader_thread(ktp_t* p, BatchQueue& queue, PipelineSync& sync)
{
    prepare_reader_processor_key();
    pipeline_thread_check(pthread_setspecific(reader_processor_key, &sync));
    for (;;) {
        queue.wait_for_space(pipeline_waits.reader_full);
        DataType item = p->func(p->shared, 0, nullptr);
        if (!item) break;
        queue.push(item);
    }
    queue.close();
    pipeline_thread_check(pthread_setspecific(reader_processor_key, nullptr));
}

static void processor_thread(ktp_t* p, BatchQueue& input, BatchQueue& output)
{
    DataType item;
    while (input.pop(item, pipeline_waits.processor_empty, true)) {
        /* Reserve before touching a SAM ring slot, not after producing it.
         * The writer's active item is outside the queue and still owns a slot. */
        output.wait_for_space(pipeline_waits.processor_full, true);
        output.push(p->func(p->shared, 1, item));
    }
    output.close();
}

static void writer_thread(ktp_t* p, BatchQueue& queue, PipelineSync& sync)
{
    prepare_reader_processor_key();
    pipeline_thread_check(pthread_setspecific(reader_processor_key, &sync));
    DataType item;
    while (queue.pop(item, pipeline_waits.writer_empty)) p->func(p->shared, 2, item);
    p->func(p->shared, 3, nullptr);
    pipeline_thread_check(pthread_setspecific(reader_processor_key, nullptr));
    sync.finish_writer();
}

extern "C" void kt_pipeline_wait_report(void)
{
    if (!pipeline_waits.enabled) return;
    const QueueWait* rows[] = {&pipeline_waits.reader_full, &pipeline_waits.processor_empty,
                              &pipeline_waits.processor_full, &pipeline_waits.writer_empty};
    const char* labels[] = {"reader waits for input queue space", "processor waits for input",
                            "processor waits for SAM buffer space", "writer waits for output"};
    fprintf(stderr, "\n  Pipeline queue waits (outside stage callback timers)\n");
    for (size_t i = 0; i < sizeof(rows) / sizeof(rows[0]); ++i)
        fprintf(stderr, "    %-43s %10.3f s  (%llu waits)\n", labels[i], rows[i]->seconds,
                (unsigned long long)rows[i]->events);
    fprintf(stderr, "    Concurrent waits overlap; do not add them to stage times.\n");
}

extern "C" void kt_pipeline_queue(int n_threads, void* (*func)(void*, int, void*), void* shared_data, int n_steps) {

    ktp_t aux;
    aux.func = func;
    aux.shared = shared_data;
    aux.n_steps = n_steps;
    assert(n_threads == 3);
    assert(n_steps == 3);

    int input_limit = read_queue_item_limit;
    if (const char *setting = std::getenv("SWBWA_INPUT_QUEUE_CHUNKS")) {
        char *end = nullptr;
        const long value = std::strtol(setting, &end, 10);
        if (!*setting || *end || value < 1 || value > 8) {
            std::fprintf(stderr, "[E::pipeline] SWBWA_INPUT_QUEUE_CHUNKS must be 1..8\n");
            std::abort();
        }
        input_limit = int(value);
    }
    std::fprintf(stderr, "[pipeline] input_queue_chunks=%d\n", input_limit);
    PipelineSync sync;
    BatchQueue read_queue(input_limit, sync), write_queue(write_queue_item_limit, sync);
    pipeline_waits = {};
    pipeline_waits.enabled = true;
    std::thread reader(reader_thread, &aux, std::ref(read_queue), std::ref(sync));
    std::thread writer(writer_thread, &aux, std::ref(write_queue), std::ref(sync));

    processor_thread(&aux, read_queue, write_queue);
    sync.drain();
    reader.join();
    writer.join();
}
