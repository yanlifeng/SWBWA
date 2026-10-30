#ifndef SWBWA_CODEC_RUNTIME_H
#define SWBWA_CODEC_RUNTIME_H

#include "swbwa_config.h"
#include <stdio.h>

#if SWBWA_ENABLE_CODEC
#ifdef __cplusplus
extern "C" {
#endif

/* Non-recursive, process-wide ownership of the CPE runtime and task state.
 * Hold across the entire alignment transaction, including host SAM layout.
 * Do not hold this lock while waiting on pipeline queues. */
void swbwa_cpe_runtime_lock(void);
void swbwa_cpe_runtime_unlock(void);

/* Cumulative alignment-only mutex wait in monotonic seconds, included in stage2.
 * Takes the runtime lock for a snapshot; call without owning it, after pipeline
 * completion for final reporting. Codec launches and this getter do not add wait. */
double swbwa_alignment_runtime_wait_seconds(void);

/* Synchronous launch of a swbwa_codec_batch from src/gzip/codec_job.h.
 * Requires completed runtime/pool initialization in main_mem. Standard modes
 * use spawn/join on the processor/main thread; reader/writer requests in all
 * modes use the shared pipeline dispatcher by default. In cgs_cross,
 * SWBWA_GZIP_CPE_DISPATCH=reader lets decoding callers publish to the persistent
 * dispatcher under the runtime lock, after initialization on the main thread.
 * With gzip input enabled, hybrid may instead decode on the available MPE helper
 * pool while CPEs are occupied; both decoding backends use independent buffers.
 * Compression still uses the processor; standard modes retain main-thread launches.
 * Do not launch codecs from arbitrary helper threads.
 * Acquires the lock internally; callers must NOT already own it.
 * The batch and every referenced buffer must remain live until return.
 * Optional timings separate dispatch/mutex wait from backend execution.
 * This function never initializes/resets the pool or halts the CPE runtime. */
void swbwa_codec_runtime_launch(void *batch, double *run_s, double *wait_s);
/* Call after all pipeline threads have joined. */
void swbwa_codec_runtime_report(FILE *stream);

#ifdef __cplusplus
}
#endif
#endif /* SWBWA_ENABLE_CODEC */

#endif /* SWBWA_CODEC_RUNTIME_H */
