#ifndef SWBWA_HOST_WORKERS_H
#define SWBWA_HOST_WORKERS_H

#include "swbwa_config.h"

typedef void (*swbwa_host_work_fn)(void *data, int worker, int workers);

#if SWBWA_HOST_MPE_THREADS > 1
int swbwa_host_workers_init(void);
void swbwa_host_workers_run(swbwa_host_work_fn function, void *data);
/* Do not wait behind I/O. A busy pool executes all logical slices on the caller. */
void swbwa_host_workers_run_ready(swbwa_host_work_fn function, void *data);
void swbwa_host_workers_destroy(void);
#else
static inline int swbwa_host_workers_init(void) { return 0; }
static inline void swbwa_host_workers_run(swbwa_host_work_fn function, void *data)
{
    function(data, 0, 1);
}
static inline void swbwa_host_workers_destroy(void) {}
static inline void swbwa_host_workers_run_ready(swbwa_host_work_fn function, void *data)
{
    function(data, 0, 1);
}
#endif

#endif
