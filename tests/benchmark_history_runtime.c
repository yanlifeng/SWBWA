/* Test-only setup equivalent to main_mem, without loading a reference index.
 * The application's backend_swbwa and cross dispatcher do all codec launches. */
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "swbwa_runtime.h"
#include "bwamem.h"
#include "swbwa_cpe.h"
#include "swbwa_host_workers.h"

extern void SLAVE_FUN(state_init());

void swbwa_benchmark_runtime_init(void)
{
    swbwa_runtime_init();
    if (swbwa_host_workers_init() != 0) MPI_Abort(MPI_COMM_WORLD, 1);
#if SWBWA_CPE_ALLOC_MODE == SWBWA_CPE_ALLOC_POOL
    swbwa_cpe_pool_params_t params;
    params.buffer = malloc(SWBWA_CPE_POOL_TOTAL_BYTES);
    if (!params.buffer) MPI_Abort(MPI_COMM_WORLD, 1);
    memset(params.buffer, 0, SWBWA_CPE_POOL_TOTAL_BYTES);
    params.bytes_per_cpe = SWBWA_CPE_POOL_BYTES_PER_CPE;
    swbwa_cpe_run((void *)slave_state_init, &params);
    /* TLS retains these addresses until process exit, just as in main_mem. */
#endif
}

void swbwa_benchmark_runtime_finish(void)
{
    swbwa_host_workers_destroy();
}
