#include <slave.h>
#include "codec_job.h"
/* Entered by the application's standard or cross entry; ownership stays outside
 * the codec. SWBWA_CODEC_LDM=2 state allocation is in rapid_state.h. */
void swbwa_codec_worker(swbwa_codec_batch *batch)
{
    unsigned i;
    for (i = (unsigned)_MYID; i < batch->count; i += batch->workers) {
        if ((unsigned)_MYID >= batch->workers) break;
        swbwa_codec_execute(&batch->jobs[i]);
    }
}
