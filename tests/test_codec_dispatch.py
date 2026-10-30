#!/usr/bin/env python3
"""Check the actual runtime dispatch code with mocked CPE launches and real queues.

This tests thread ownership and mutual exclusion, not Sunway MMIO or athread.
"""
from pathlib import Path
import os
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

PREFIX = r'''
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "codec_job.h"
#include "swbwa_pipeline.h"
#define SWBWA_CPE_PROGRESS_CODEC 0
static pthread_t owner, reader;
static pthread_mutex_t runtime_mutex = PTHREAD_MUTEX_INITIALIZER;
#define swbwa_cpe_runtime_mutex runtime_mutex
static int bwa_verbose = 4, init_calls, decoded, encoded, written;
static int active, direct_decode, queued_decode, next;
static int mock_mpe_calls, test_hybrid, compute_started;
typedef struct { void *codec_batch; } swbwa_cpe_task_t;
typedef struct { int initialized; } cross_runtime_t;
static cross_runtime_t swbwa_shared_cross_runtime;
static swbwa_cpe_task_t parameters;
static double GetTime(void) {
    struct timespec t;
    assert(!clock_gettime(CLOCK_MONOTONIC, &t));
    return t.tv_sec + t.tv_nsec * 1e-9;
}
static void err_fatal(const char *where, const char *fmt, ...) {
    va_list ap;
    fprintf(stderr, "%s: ", where);
    va_start(ap, fmt); vfprintf(stderr, fmt, ap); va_end(ap);
    exit(42);
}
static int swbwa_mpi_rank(void) { return 0; }
static void swbwa_cpe_runtime_lock(void) { assert(!pthread_mutex_lock(&runtime_mutex)); }
static void swbwa_cpe_runtime_unlock(void) { assert(!pthread_mutex_unlock(&runtime_mutex)); }
static swbwa_cpe_task_t *swbwa_shared_params_locked(void) { return &parameters; }
static void swbwa_cross_runtime_init(cross_runtime_t *r, swbwa_cpe_task_t *p) {
    (void)p;
    if (r->initialized) return;
    assert(pthread_equal(pthread_self(), owner));
    ++init_calls; r->initialized = 1;
}
static void slave_swbwa_codec_cross(void) {}
static void slave_swbwa_codec(void) {}
static long swbwa_cross_entry(cross_runtime_t *r, void *entry) {
    assert(r->initialized); return (long)entry;
}
static void execute(void) {
    swbwa_codec_batch *batch = parameters.codec_batch;
    assert(!__atomic_exchange_n(&active, 1, __ATOMIC_SEQ_CST));
    assert(batch && batch->count == 1);
    if (batch->jobs[0].operation == SWBWA_CODEC_GZIP ||
        batch->jobs[0].operation == SWBWA_CODEC_DEFLATE) {
        assert(pthread_equal(pthread_self(), owner)); ++encoded;
    } else {
        ++decoded;
        if (pthread_equal(pthread_self(), reader)) ++direct_decode;
        else { assert(pthread_equal(pthread_self(), owner)); ++queued_decode; }
    }
    usleep(100);
    assert(__atomic_exchange_n(&active, 0, __ATOMIC_SEQ_CST));
}
static void swbwa_cross_run(long entry) { (void)entry; execute(); }
static void swbwa_cpe_run_with_mpi_progress(void *entry, swbwa_cpe_task_t *p, int phase) {
    (void)entry; (void)p; (void)phase;
    assert(pthread_equal(pthread_self(), owner)); execute();
}
int swbwa_codec_mpe_decode(swbwa_codec_batch *batch) {
    assert(pthread_equal(pthread_self(), reader));
    assert(batch->count == 1 && (batch->jobs[0].operation == SWBWA_CODEC_DISCOVER ||
                               batch->jobs[0].operation == SWBWA_CODEC_RESOLVE));
    ++decoded;
    __atomic_add_fetch(&mock_mpe_calls, 1, __ATOMIC_SEQ_CST);
    usleep(100);
    return 1;
}
'''

SUFFIX = r'''
extern void kt_pipeline_queue(int, void *(*)(void *, int, void *), void *, int);
static void launch(int operation) {
    swbwa_codec_job job = {0};
    swbwa_codec_batch batch = {0};
    double run = -1, wait = -1;
    job.operation = operation;
    batch.jobs = &job; batch.count = 1; batch.workers = SWBWA_CPE_COUNT;
    swbwa_codec_runtime_launch(&batch, &run, &wait);
    assert(run >= 0 && wait >= 0);
}
static void *process(void *shared, int step, void *item) {
    (void)shared;
    if (step == 0) {
        if (next == 0) reader = pthread_self();
        if (next == 24) return NULL;
        if (test_hybrid && next == 1)
            while (!__atomic_load_n(&compute_started, __ATOMIC_SEQ_CST)) usleep(20);
        launch(SWBWA_CODEC_DISCOVER);
        launch(SWBWA_CODEC_RESOLVE);
        launch(SWBWA_CODEC_DEFLATE);
        return (void *)(intptr_t)++next;
    }
    if (step == 1) {
        assert(pthread_equal(pthread_self(), owner));
        swbwa_cpe_runtime_lock();
        swbwa_cross_runtime_init(&swbwa_shared_cross_runtime, &parameters);
        assert(!__atomic_exchange_n(&active, 1, __ATOMIC_SEQ_CST));
        if (test_hybrid && (intptr_t)item == 1) {
            __atomic_store_n(&compute_started, 1, __ATOMIC_SEQ_CST);
            while (!__atomic_load_n(&mock_mpe_calls, __ATOMIC_SEQ_CST)) usleep(20);
        }
        usleep(500);
        assert(__atomic_exchange_n(&active, 0, __ATOMIC_SEQ_CST));
        swbwa_cpe_runtime_unlock();
        launch(SWBWA_CODEC_GZIP);
        return item;
    }
    if (step == 2) {
        assert((intptr_t)item == ++written);
        usleep(700);
    }
    return NULL;
}
int main(void) {
    owner = pthread_self();
    const char *mode = getenv("SWBWA_GZIP_CPE_DISPATCH");
    test_hybrid = mode && !strcmp(mode, "hybrid");
    kt_pipeline_queue(3, process, NULL, 3);
    assert(!parameters.codec_batch);
    assert(init_calls == 1 && decoded == 48 && encoded == 48 && written == 24);
    int direct = SWBWA_USE_CROSS_SEGMENT && mode && strcmp(mode, "processor");
    assert(direct_decode == (direct ? 48 - mock_mpe_calls : 0));
    assert(queued_decode == (direct ? 0 : 48 - mock_mpe_calls));
    assert((mock_mpe_calls > 0) == test_hybrid);
    assert(swbwa_codec_dispatch_stats[1].calls == 48);
    assert(swbwa_codec_dispatch_stats[0].calls == 48);
    assert(swbwa_codec_dispatch_stats[0].direct_calls == 0);
    assert(swbwa_codec_dispatch_stats[1].mpe_calls == (unsigned)mock_mpe_calls);
    swbwa_codec_runtime_report(stderr);
    puts("CODEC_DISPATCH_PASS");
    return 0;
}
'''


def main():
    source = (ROOT / 'src/host/bwamem.c').read_text()
    begin = source.index('static pthread_once_t swbwa_codec_dispatch_once')
    end = source.index('\n#endif\n\nvoid mem_process_seqs_merge2', begin)
    with tempfile.TemporaryDirectory(prefix='swbwa-dispatch-') as tmp:
        tmp = Path(tmp)
        harness = tmp / 'dispatch.c'
        harness.write_text(PREFIX + source[begin:end] + SUFFIX)
        for cross in (False, True):
            flags = ['-O1', '-g', '-pthread', '-fsanitize=address,undefined',
                     '-fno-sanitize-recover=all', '-fno-pie', '-I' + str(ROOT / 'include'),
                     '-I' + str(ROOT / 'src/gzip'), '-include', str(ROOT / 'include/swbwa_config.h'),
                     '-DSWBWA_ENABLE_GZIP_INPUT=1',
                     '-DSWBWA_EXEC_MODE=' + ('SWBWA_EXEC_CGS_CROSS' if cross else 'SWBWA_EXEC_SINGLE_CG')]
            obj = tmp / 'dispatch.o'
            subprocess.run([os.getenv('CC', 'cc'), '-std=gnu11', *flags,
                            '-c', str(harness), '-o', str(obj)], check=True)
            binary = tmp / 'dispatch'
            subprocess.run([os.getenv('CXX', 'c++'), '-std=c++11', *flags, '-no-pie',
                            str(obj), str(ROOT / 'src/host/kthread.cpp'), '-o', str(binary)], check=True)
            for mode in (None, 'processor', 'reader', 'hybrid', 'bad'):
                env = dict(os.environ, ASAN_OPTIONS='detect_leaks=0:halt_on_error=1')
                env.pop('SWBWA_GZIP_CPE_DISPATCH', None)
                if mode is not None:
                    env['SWBWA_GZIP_CPE_DISPATCH'] = mode
                result = subprocess.run([str(binary)], env=env, capture_output=True, text=True, timeout=15)
                assert result.returncode == (42 if mode == 'bad' else 0), result.stderr
                if mode != 'bad':
                    assert 'CODEC_DISPATCH_PASS' in result.stdout, result.stderr
                print('PASS', 'cross' if cross else 'single', mode)


if __name__ == '__main__':
    main()
