/* Link with swbwa_mpi.c using -ffunction-sections -Wl,--gc-sections.
 * Each fatal case MUST fail promptly as a whole MPI job, never time out.
 * "inspect" only reports the native policy and exits normally; it does not
 * intentionally enter the vendor's debug-hang path. No CPEs are launched. */
#include "swbwa_mpi.h"
#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef __sw_host__
extern int swch_hang_when_error __attribute__((weak));
#endif

static MPI_Comm fault_comm;

static void *fail(void *unused)
{
    (void)unused;
    fprintf(stderr, "FATAL_EXIT_TEST rank=%d entering MPI_Abort\n", swbwa_mpi_rank());
    fflush(stderr);
    MPI_Abort(fault_comm, 7);
    _exit(7);
}

int main(int argc, char **argv)
{
    int provided = 0, rank = 0, ranks = 0;
    const char *mode = argc == 2 ? argv[1] : "";
    int inspect = strcmp(mode, "inspect") == 0;
    int all = strcmp(mode, "all") == 0;
    int threaded = strcmp(mode, "thread") == 0;
    int normal = strcmp(mode, "normal") == 0;
    if (!inspect && !all && !threaded && !normal &&
        strcmp(mode, "world") && strcmp(mode, "dup")) {
        fprintf(stderr, "usage: %s inspect|normal|world|dup|all|thread\n", argv[0]);
        return 2;
    }
    if (inspect) {
        if (MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided) != MPI_SUCCESS)
            return 2;
    } else if (swbwa_mpi_init(&argc, &argv) != 0) {
        return 2;
    }
#ifdef __sw_host__
    if (!inspect && (!getenv("SWCH_HANG_WHEN_ERROR") ||
                    strcmp(getenv("SWCH_HANG_WHEN_ERROR"), "0")))
        swbwa_mpi_abort("SWCH error-hang policy was not disabled");
#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Query_thread(&provided);
    fprintf(stderr, "FATAL_EXIT_TEST rank=%d/%d mode=%s provided=%d", rank, ranks, mode, provided);
#ifdef __sw_host__
    if (&swch_hang_when_error) fprintf(stderr, " native_hang=%d", swch_hang_when_error);
#endif
    fprintf(stderr, "\n");
    fflush(stderr);
    if (inspect || normal) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) puts("FATAL_EXIT_TEST normal shutdown");
        if (inspect) MPI_Finalize();
        else swbwa_mpi_finalize();
        return 0;
    }
    fault_comm = MPI_COMM_WORLD;
    if (strcmp(mode, "world") && MPI_Comm_dup(MPI_COMM_WORLD, &fault_comm) != MPI_SUCCESS)
        swbwa_mpi_abort("cannot duplicate test communicator");
    MPI_Barrier(MPI_COMM_WORLD);
    if (all || rank == ranks - 1) {
        if (threaded) {
            pthread_t thread;
            if (pthread_create(&thread, NULL, fail, NULL) != 0)
                swbwa_mpi_abort("cannot create test thread");
            pthread_join(thread, NULL);
        } else {
            fail(NULL);
        }
    }
    /* Other ranks wait in a different communicator, as in a failed Reader. */
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "FATAL_EXIT_TEST unexpected return from failed job\n");
    _exit(3);
}
