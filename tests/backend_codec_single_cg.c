/* Standalone benchmark backend: one native 64-CPE group per MPI rank. */
#include <athread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "codec_job.h"

extern void slave_swbwa_codec_worker(void *);

void swbwa_codec_init(void) { athread_init(); }
void swbwa_codec_finish(void) {}

void swbwa_codec_launch(swbwa_codec_batch *batch)
{
    struct timespec begin, end;
    if (!batch->workers || batch->workers > 64) {
        fprintf(stderr, "single-CG codec requires 1..64 workers\n");
        abort();
    }
    clock_gettime(CLOCK_MONOTONIC, &begin);
    __real_athread_spawn((void *)slave_swbwa_codec_worker, batch, 1);
    athread_join();
    clock_gettime(CLOCK_MONOTONIC, &end);
    batch->run_s = end.tv_sec - begin.tv_sec + (end.tv_nsec - begin.tv_nsec) * 1e-9;
    batch->wait_s = 0;
}
