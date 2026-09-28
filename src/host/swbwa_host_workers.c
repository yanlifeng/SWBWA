#include "swbwa_host_workers.h"

#if SWBWA_HOST_MPE_THREADS > 1
#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#ifdef __sw_host__
#include <athread.h>
#endif

typedef struct {
    pthread_t thread;
    int worker, target, actual, bind_result;
} host_worker_t;

static struct {
    pthread_mutex_t dispatch;
    pthread_mutex_t mutex;
    pthread_cond_t work, done;
    host_worker_t workers[SWBWA_HOST_MPE_THREADS - 1];
    swbwa_host_work_fn function;
    void *data;
    unsigned long generation;
    int created, ready, remaining, stopping, initialized;
} pool = {
    .dispatch = PTHREAD_MUTEX_INITIALIZER,
    .mutex = PTHREAD_MUTEX_INITIALIZER,
    .work = PTHREAD_COND_INITIALIZER,
    .done = PTHREAD_COND_INITIALIZER
};

static int current_cg(void)
{
#ifdef __sw_host__
    unsigned long id;
    __asm__ volatile("rcid %0" : "=r"(id));
    return (int)(id & 7);
#else
    return -1; /* Native tests exercise synchronization, not Sunway placement. */
#endif
}

static void *worker_entry(void *opaque)
{
    host_worker_t *worker = opaque;
    unsigned long generation = 0;
#ifdef __sw_host__
    worker->bind_result = SET_TO_SWCG(worker->target);
#endif
    worker->actual = current_cg();
    pthread_mutex_lock(&pool.mutex);
    ++pool.ready;
    pthread_cond_broadcast(&pool.done);
    while (!pool.stopping) {
        swbwa_host_work_fn function;
        void *data;
        while (!pool.stopping && generation == pool.generation)
            pthread_cond_wait(&pool.work, &pool.mutex);
        if (pool.stopping) break;
        generation = pool.generation;
        function = pool.function;
        data = pool.data;
        pthread_mutex_unlock(&pool.mutex);
        function(data, worker->worker, SWBWA_HOST_MPE_THREADS);
        pthread_mutex_lock(&pool.mutex);
        if (--pool.remaining == 0) pthread_cond_signal(&pool.done);
    }
    pthread_mutex_unlock(&pool.mutex);
    return NULL;
}

int swbwa_host_workers_init(void)
{
    int i, error = 0, home = current_cg();
    if (pool.initialized) return 0;
#ifdef __sw_host__
    if (home < 0 || home >= SWBWA_HOST_MPE_THREADS) {
        errno = EINVAL;
        return -1;
    }
#endif
    for (i = 1; i < SWBWA_HOST_MPE_THREADS; ++i) {
        host_worker_t *worker = &pool.workers[i - 1];
        worker->worker = i;
        worker->target = (home + i) % SWBWA_HOST_MPE_THREADS;
        error = pthread_create(&worker->thread, NULL, worker_entry, worker);
        if (error != 0) break;
        ++pool.created;
    }
    pthread_mutex_lock(&pool.mutex);
    while (pool.ready < pool.created) pthread_cond_wait(&pool.done, &pool.mutex);
    pthread_mutex_unlock(&pool.mutex);
    for (i = 0; i < pool.created; ++i) {
        host_worker_t *worker = &pool.workers[i];
#ifdef __sw_host__
        if (worker->actual != worker->target) error = EINVAL;
#endif
        fprintf(stderr, "[MPE helper] worker=%d target_cg=%d actual_cg=%d bind_result=%d\n",
                worker->worker, worker->target, worker->actual, worker->bind_result);
    }
    if (error != 0) {
        swbwa_host_workers_destroy();
        errno = error;
        return -1;
    }
    pool.initialized = 1;
    fprintf(stderr, "[MPE helpers] threads=%d caller_cg=%d placement=%s\n",
            SWBWA_HOST_MPE_THREADS, home, home < 0 ? "native-unbound" : "verified");
    return 0;
}

static void run_locked(swbwa_host_work_fn function, void *data)
{
    if (!pool.initialized) {
        function(data, 0, 1);
        pthread_mutex_unlock(&pool.dispatch);
        return;
    }
    pthread_mutex_lock(&pool.mutex);
    pool.function = function;
    pool.data = data;
    pool.remaining = pool.created;
    ++pool.generation;
    pthread_cond_broadcast(&pool.work);
    pthread_mutex_unlock(&pool.mutex);
    function(data, 0, SWBWA_HOST_MPE_THREADS);
    pthread_mutex_lock(&pool.mutex);
    while (pool.remaining != 0) pthread_cond_wait(&pool.done, &pool.mutex);
    pool.function = NULL;
    pool.data = NULL;
    pthread_mutex_unlock(&pool.mutex);
    pthread_mutex_unlock(&pool.dispatch);
}

void swbwa_host_workers_run(swbwa_host_work_fn function, void *data)
{
    pthread_mutex_lock(&pool.dispatch);
    run_locked(function, data);
}

void swbwa_host_workers_run_ready(swbwa_host_work_fn function, void *data)
{
    if (pthread_mutex_trylock(&pool.dispatch) == 0) {
        run_locked(function, data);
    } else {
        int worker;
        /* Keep slice geometry identical across measure/pack calls even when
         * just one of them obtains the pool. Never borrow its task storage. */
        for (worker = 0; worker < SWBWA_HOST_MPE_THREADS; ++worker)
            function(data, worker, SWBWA_HOST_MPE_THREADS);
    }
}

void swbwa_host_workers_destroy(void)
{
    int i;
    pthread_mutex_lock(&pool.dispatch);
    pthread_mutex_lock(&pool.mutex);
    pool.stopping = 1;
    pthread_cond_broadcast(&pool.work);
    pthread_mutex_unlock(&pool.mutex);
    for (i = 0; i < pool.created; ++i) pthread_join(pool.workers[i].thread, NULL);
    pool.created = pool.ready = pool.remaining = 0;
    pool.stopping = pool.initialized = 0;
    pool.generation = 0;
    pthread_mutex_unlock(&pool.dispatch);
}
#endif
