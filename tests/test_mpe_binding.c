/* Sunway host-only probe. Run on one exclusively allocated six-CG node. */
#include <athread.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { WORKERS = 6, WORDS = 1 << 20, ROUNDS = 100 };
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condition = PTHREAD_COND_INITIALIZER;
static int ready, release_workers;
static __thread int tls_value;
typedef struct {
    int target, actual, bind_result, ok;
    uint64_t sum;
    uint64_t *data;
} worker_t;
static worker_t workers[WORKERS];

static int current_cg(void)
{
    unsigned long id;
    __asm__ volatile("rcid %0" : "=r"(id));
    return (int)(id & 7);
}

static void work(worker_t *w)
{
    int round, i;
    volatile unsigned long stack_value = (unsigned long)w->target + 100;
    tls_value = w->target + 1000;
    for (round = 0; round < ROUNDS; ++round) {
        uint64_t sum = 0;
        for (i = 0; i < WORDS; ++i) {
            w->data[i] = (uint64_t)i + (unsigned)w->target + (unsigned)round;
            sum += w->data[i];
        }
        w->sum = sum;
    }
    w->ok = tls_value == w->target + 1000 &&
            stack_value == (unsigned long)w->target + 100 &&
            current_cg() == w->target;
}

static void *entry(void *opaque)
{
    worker_t *w = opaque;
    w->bind_result = SET_TO_SWCG(w->target);
    w->actual = current_cg();
    pthread_mutex_lock(&mutex);
    ++ready;
    pthread_cond_broadcast(&condition);
    while (!release_workers) pthread_cond_wait(&condition, &mutex);
    pthread_mutex_unlock(&mutex);
    if (w->actual == w->target) work(w);
    return NULL;
}

int main(void)
{
    pthread_t threads[WORKERS - 1];
    int home = current_cg(), i, created = 0, failed = 0;
    if (home < 0 || home >= WORKERS) return 1;
    printf("primary CG=%d\n", home);
    fflush(stdout);
    for (i = 0; i < WORKERS; ++i) {
        workers[i].target = (home + i) % WORKERS;
        workers[i].data = malloc((size_t)WORDS * sizeof(uint64_t));
        if (workers[i].data == NULL) return 2;
    }
    workers[0].actual = home;
    for (i = 1; i < WORKERS; ++i) {
        int error = pthread_create(&threads[i - 1], NULL, entry, &workers[i]);
        if (error) { fprintf(stderr, "pthread_create error=%d\n", error); failed = 1; break; }
        ++created;
    }
    pthread_mutex_lock(&mutex);
    while (ready < created) pthread_cond_wait(&condition, &mutex);
    for (i = 1; i <= created; ++i)
        printf("helper=%d target=%d actual=%d bind_result=%d\n", i,
               workers[i].target, workers[i].actual, workers[i].bind_result);
    fflush(stdout);
    release_workers = 1;
    pthread_cond_broadcast(&condition);
    pthread_mutex_unlock(&mutex);
    work(&workers[0]);
    for (i = 0; i < created; ++i) pthread_join(threads[i], NULL);
    for (i = 0; i < WORKERS; ++i) {
        uint64_t expected = (uint64_t)WORDS * (WORDS - 1) / 2 +
                            (uint64_t)WORDS * (workers[i].target + ROUNDS - 1);
        int ok = workers[i].ok && workers[i].sum == expected;
        printf("worker=%d CG=%d checksum=%llu expected=%llu ok=%d\n", i,
               workers[i].actual, (unsigned long long)workers[i].sum,
               (unsigned long long)expected, ok);
        if (!ok) failed = 1;
        free(workers[i].data);
    }
    printf("MPE binding probe: %s\n", failed ? "FAIL" : "PASS");
    return failed;
}
