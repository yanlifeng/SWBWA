#ifndef SWBWA_PIPELINE_H
#define SWBWA_PIPELINE_H

#ifdef __cplusplus
extern "C" {
#endif

/* A reader callback can synchronously request work on the pipeline's processor
 * (calling/main) thread. Other callers execute directly. The argument remains
 * owned by the caller and must live until return. The callback must not wait on
 * pipeline queues or recursively request processor work. */
void swbwa_pipeline_call_on_processor(void (*call)(void *), void *argument);

#ifdef __cplusplus
}
#endif

#endif
