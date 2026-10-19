#ifndef SWBWA_CODEC_RUNTIME_H
#define SWBWA_CODEC_RUNTIME_H

#include "swbwa_config.h"

#if SWBWA_ENABLE_GZIP_INPUT
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
 * use spawn/join on the processor/main thread; reader requests are handed off
 * through the input queue. cgs_cross uses the relocated persistent dispatcher.
 * Call from the main thread or pipeline reader, not an arbitrary helper thread.
 * Acquires the lock internally; callers must NOT already own it.
 * The batch and every referenced buffer must remain live until return.
 * Optional timings separate dispatch/mutex wait from runtime initialization/execution.
 * This function never initializes/resets the pool or halts the CPE runtime. */
void swbwa_codec_runtime_launch(void *batch, double *run_s, double *wait_s);

#ifdef __cplusplus
}
#endif
#endif /* SWBWA_ENABLE_GZIP_INPUT */

#endif /* SWBWA_CODEC_RUNTIME_H */
