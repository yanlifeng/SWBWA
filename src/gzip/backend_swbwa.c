#include "codec_job.h"
#include "swbwa_codec_runtime.h"

/* main_mem owns initialization, pools and lifetime of the shared runtime. */
void swbwa_codec_init(void) {}
void swbwa_codec_finish(void) {}

void swbwa_codec_launch(swbwa_codec_batch *batch)
{
    swbwa_codec_runtime_launch(batch, &batch->run_s, &batch->wait_s);
}
