#include <assert.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include "swbwa_host_workers.h"

enum { COUNT = 1001, ROUNDS = 200 };
static pthread_mutex_t held_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t held_changed = PTHREAD_COND_INITIALIZER;
static int entered, released;
typedef struct { unsigned values[COUNT], seed; } task_t;
static void fill(void *opaque, int worker, int workers)
{
    task_t *task = opaque;
    size_t i, start = (size_t)COUNT * worker / workers;
    size_t end = (size_t)COUNT * (worker + 1) / workers;
    for (i = start; i < end; ++i) task->values[i] = (unsigned)i + task->seed;
}
static void *caller(void *opaque)
{
    task_t task = {{0}, (unsigned)(uintptr_t)opaque};
    int round, i;
    for (round = 0; round < ROUNDS; ++round) {
        if (round & 1) swbwa_host_workers_run_ready(fill, &task);
        else swbwa_host_workers_run(fill, &task);
        for (i = 0; i < COUNT; ++i) assert(task.values[i] == (unsigned)i + task.seed);
        ++task.seed;
    }
    return NULL;
}

static void hold_pool(void *opaque, int worker, int workers)
{
    (void)opaque;
    (void)workers;
    if (worker != 0) return;
    pthread_mutex_lock(&held_mutex);
    entered = 1;
    pthread_cond_broadcast(&held_changed);
    while (!released) pthread_cond_wait(&held_changed, &held_mutex);
    pthread_mutex_unlock(&held_mutex);
}

static void *blocked_reader(void *opaque)
{
    swbwa_host_workers_run(hold_pool, opaque);
    return NULL;
}

static void check_busy_fallback(void)
{
    pthread_t reader;
    task_t task = {{0}, 9000};
    int i;
    entered = released = 0;
    assert(pthread_create(&reader, NULL, blocked_reader, NULL) == 0);
    pthread_mutex_lock(&held_mutex);
    while (!entered) pthread_cond_wait(&held_changed, &held_mutex);
    pthread_mutex_unlock(&held_mutex);
    swbwa_host_workers_run_ready(fill, &task);
    for (i = 0; i < COUNT; ++i) assert(task.values[i] == (unsigned)i + task.seed);
    pthread_mutex_lock(&held_mutex);
    released = 1;
    pthread_cond_broadcast(&held_changed);
    pthread_mutex_unlock(&held_mutex);
    assert(pthread_join(reader, NULL) == 0);
}
int main(void)
{
    int cycle;
    for (cycle = 0; cycle < 3; ++cycle) {
        pthread_t other[2];
        assert(swbwa_host_workers_init() == 0);
        check_busy_fallback();
        assert(pthread_create(&other[0], NULL, caller, (void *)(uintptr_t)1234) == 0);
        assert(pthread_create(&other[1], NULL, caller, (void *)(uintptr_t)5678) == 0);
        caller(NULL);
        assert(pthread_join(other[0], NULL) == 0);
        assert(pthread_join(other[1], NULL) == 0);
        swbwa_host_workers_destroy();
    }
    puts("host worker concurrent dispatch and restart: PASS");
    return 0;
}
