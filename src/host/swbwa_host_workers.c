#include "swbwa_host_workers.h"

#if SWBWA_HOST_MPE_THREADS > 1
#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#ifdef __sw_host__
#include <athread.h>
#endif

typedef struct host_worker_pool host_worker_pool_t;
typedef struct {
    pthread_t thread;
    int worker, target, actual, bind_result;
    host_worker_pool_t *pool;
} host_worker_t;

struct host_worker_pool {
    pthread_mutex_t dispatch;
    pthread_mutex_t mutex;
    pthread_cond_t work, done;
    host_worker_t workers[SWBWA_HOST_MPE_THREADS - 1];
    swbwa_host_work_fn function;
    void *data;
    unsigned long generation;
    int created, ready, remaining, stopping, initialized;
    const char *name;
};

static host_worker_pool_t compute_pool = {
    .dispatch = PTHREAD_MUTEX_INITIALIZER,
    .mutex = PTHREAD_MUTEX_INITIALIZER,
    .work = PTHREAD_COND_INITIALIZER,
    .done = PTHREAD_COND_INITIALIZER,
    .name = "MPE helpers"
};
static host_worker_pool_t io_pool = {
    .dispatch = PTHREAD_MUTEX_INITIALIZER,
    .mutex = PTHREAD_MUTEX_INITIALIZER,
    .work = PTHREAD_COND_INITIALIZER,
    .done = PTHREAD_COND_INITIALIZER,
    .name = "MPE I/O helpers"
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
    host_worker_pool_t *pool = worker->pool;
    unsigned long generation = 0;
#ifdef __sw_host__
    worker->bind_result = SET_TO_SWCG(worker->target);
#endif
    worker->actual = current_cg();
    pthread_mutex_lock(&pool->mutex);
    ++pool->ready;
    pthread_cond_broadcast(&pool->done);
    while (!pool->stopping) {
        swbwa_host_work_fn function;
        void *data;
        while (!pool->stopping && generation == pool->generation)
            pthread_cond_wait(&pool->work, &pool->mutex);
        if (pool->stopping) break;
        generation = pool->generation;
        function = pool->function;
        data = pool->data;
        pthread_mutex_unlock(&pool->mutex);
        function(data, worker->worker, SWBWA_HOST_MPE_THREADS);
        pthread_mutex_lock(&pool->mutex);
        if (--pool->remaining == 0) pthread_cond_signal(&pool->done);
    }
    pthread_mutex_unlock(&pool->mutex);
    return NULL;
}

/* Called with this pool's dispatch lock held and no callback in flight. */
static void destroy_locked(host_worker_pool_t *pool)
{
    int i;
    pthread_mutex_lock(&pool->mutex);
    pool->stopping = 1;
    pthread_cond_broadcast(&pool->work);
    pthread_mutex_unlock(&pool->mutex);
    for (i = 0; i < pool->created; ++i) pthread_join(pool->workers[i].thread, NULL);
    pool->created = pool->ready = pool->remaining = 0;
    pool->stopping = pool->initialized = 0;
    pool->generation = 0;
}

static int init_locked(host_worker_pool_t *pool)
{
    int i, error = 0, home = current_cg();
    if (pool->initialized) return 0;
#ifdef __sw_host__
    if (home < 0 || home >= SWBWA_HOST_MPE_THREADS) {
        errno = EINVAL;
        return -1;
    }
#endif
    for (i = 1; i < SWBWA_HOST_MPE_THREADS; ++i) {
        host_worker_t *worker = &pool->workers[i - 1];
        worker->pool = pool;
        worker->worker = i;
        worker->target = (home + i) % SWBWA_HOST_MPE_THREADS;
        error = pthread_create(&worker->thread, NULL, worker_entry, worker);
        if (error != 0) break;
        ++pool->created;
    }
    pthread_mutex_lock(&pool->mutex);
    while (pool->ready < pool->created) pthread_cond_wait(&pool->done, &pool->mutex);
    pthread_mutex_unlock(&pool->mutex);
    for (i = 0; i < pool->created; ++i) {
        host_worker_t *worker = &pool->workers[i];
#ifdef __sw_host__
        if (worker->actual != worker->target) error = EINVAL;
#endif
        fprintf(stderr, "[%s] worker=%d target_cg=%d actual_cg=%d bind_result=%d\n",
                pool->name, worker->worker, worker->target, worker->actual, worker->bind_result);
    }
    if (error != 0) {
        destroy_locked(pool);
        errno = error;
        return -1;
    }
    pool->initialized = 1;
    fprintf(stderr, "[%s] threads=%d caller_cg=%d placement=%s\n",
            pool->name, SWBWA_HOST_MPE_THREADS, home, home < 0 ? "native-unbound" : "verified");
    return 0;
}

int swbwa_host_workers_init(void)
{
    int status;
    pthread_mutex_lock(&compute_pool.dispatch);
    status = init_locked(&compute_pool);
    pthread_mutex_unlock(&compute_pool.dispatch);
    return status;
}

static void run_locked(host_worker_pool_t *pool, swbwa_host_work_fn function, void *data)
{
    if (!pool->initialized) {
        function(data, 0, 1);
        pthread_mutex_unlock(&pool->dispatch);
        return;
    }
    pthread_mutex_lock(&pool->mutex);
    pool->function = function;
    pool->data = data;
    pool->remaining = pool->created;
    ++pool->generation;
    pthread_cond_broadcast(&pool->work);
    pthread_mutex_unlock(&pool->mutex);
    function(data, 0, SWBWA_HOST_MPE_THREADS);
    pthread_mutex_lock(&pool->mutex);
    while (pool->remaining != 0) pthread_cond_wait(&pool->done, &pool->mutex);
    pool->function = NULL;
    pool->data = NULL;
    pthread_mutex_unlock(&pool->mutex);
    pthread_mutex_unlock(&pool->dispatch);
}

void swbwa_host_workers_run(swbwa_host_work_fn function, void *data)
{
    pthread_mutex_lock(&compute_pool.dispatch);
    run_locked(&compute_pool, function, data);
}

void swbwa_host_workers_run_ready(swbwa_host_work_fn function, void *data)
{
    if (pthread_mutex_trylock(&compute_pool.dispatch) == 0) {
        run_locked(&compute_pool, function, data);
    } else {
        int worker;
        /* Keep slice geometry identical across measure/pack calls even when
         * just one of them obtains the pool. Never borrow its task storage. */
        for (worker = 0; worker < SWBWA_HOST_MPE_THREADS; ++worker)
            function(data, worker, SWBWA_HOST_MPE_THREADS);
    }
}

int swbwa_host_workers_try_run(swbwa_host_work_fn function, void *data)
{
    if (pthread_mutex_trylock(&compute_pool.dispatch) != 0) return 0;
    run_locked(&compute_pool, function, data);
    return 1;
}

int swbwa_host_io_workers_run(swbwa_host_work_fn function, void *data)
{
    pthread_mutex_lock(&io_pool.dispatch);
    if (init_locked(&io_pool) != 0) {
        pthread_mutex_unlock(&io_pool.dispatch);
        return -1;
    }
    run_locked(&io_pool, function, data);
    return 0;
}

void swbwa_host_workers_destroy(void)
{
    pthread_mutex_lock(&io_pool.dispatch);
    destroy_locked(&io_pool);
    pthread_mutex_unlock(&io_pool.dispatch);
    pthread_mutex_lock(&compute_pool.dispatch);
    destroy_locked(&compute_pool);
    pthread_mutex_unlock(&compute_pool.dispatch);
}
#endif
