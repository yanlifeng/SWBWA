#ifndef SWBWA_RUNTIME_H
#define SWBWA_RUNTIME_H

#include <athread.h>

#include "swbwa_config.h"

static inline void swbwa_runtime_init(void)
{
#if SWBWA_USE_CGS
    athread_init_cgs();
#else
    athread_init();
#endif
}

static inline void swbwa_cpe_spawn(void *entry, void *argument)
{
#if SWBWA_USE_CGS
    __real_athread_spawn_cgs(entry, argument, 1);
#else
    __real_athread_spawn(entry, argument, 1);
#endif
}

static inline void swbwa_cpe_join(void)
{
#if SWBWA_USE_CGS
    athread_join_cgs();
#else
    athread_join();
#endif
}

static inline void swbwa_cpe_run(void *entry, void *argument)
{
    swbwa_cpe_spawn(entry, argument);
    swbwa_cpe_join();
}

#endif /* SWBWA_RUNTIME_H */
