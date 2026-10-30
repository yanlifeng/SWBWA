#ifndef SWBWA_PIPELINE_H
#define SWBWA_PIPELINE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Reader/writer callbacks can synchronously request work on the processor
 * (calling/main) thread. Other callers execute directly. The argument remains
 * owned by the caller and must live until return. The callback must not wait on
 * pipeline queues or recursively request processor work. */
void swbwa_pipeline_call_on_processor(void (*call)(void *), void *argument);

/* A nested producer may inherit the caller's handoff context. It must clear
 * that context and join before its parent pipeline callback reaches EOF. */
void *swbwa_pipeline_context(void);
void swbwa_pipeline_set_context(void *context);

#ifdef __cplusplus
}
#endif

#endif
