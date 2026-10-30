#include "codec_job.h"
#include "swbwa_host_workers.h"

void swbwa_codec_execute_mpe(swbwa_codec_job *job);

typedef struct {
    swbwa_codec_batch *batch;
    unsigned next;
} mpe_decode_task;

static void decode_jobs(void *opaque, int worker, int workers)
{
    mpe_decode_task *task = opaque;
    unsigned index;
    (void)worker;
    (void)workers;
    while ((index = __sync_fetch_and_add(&task->next, 1)) < task->batch->count)
        swbwa_codec_execute_mpe(&task->batch->jobs[index]);
}

int swbwa_codec_mpe_decode(swbwa_codec_batch *batch)
{
    mpe_decode_task task = {batch, 0};
    unsigned i;
    for (i = 0; i < batch->count; ++i)
        if (batch->jobs[i].operation != SWBWA_CODEC_DISCOVER &&
            batch->jobs[i].operation != SWBWA_CODEC_RESOLVE) return 0;
    return swbwa_host_workers_try_run(decode_jobs, &task);
}
