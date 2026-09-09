#define _POSIX_C_SOURCE 200809L

#include <stdint.h>

struct cpe_progress_task {
    volatile uint64_t done;
    uint64_t iterations;
    volatile uint64_t checksum;
};

#ifdef TEST_MPI_RMA_BUILD_CPE_KERNEL

#include <slave.h>

void test_mpi_rma_cpe_work(struct cpe_progress_task *task)
{
    uint64_t value = UINT64_C(0x9e3779b97f4a7c15) + (uint64_t)_MYID;
    uint64_t i;

    for (i = 0; i < task->iterations; ++i)
        value = value * UINT64_C(2862933555777941757) +
                UINT64_C(3037000493);

    athread_ssync_array();
    if (_MYID == 0) {
        task->checksum = value;
        asm volatile("memb\n\t" ::: "memory");
        task->done = 1;
        flush_slave_cache();
    }
}

#else

#include <mpi.h>

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef TEST_MPI_REQUIRED_THREAD_LEVEL
#define TEST_MPI_REQUIRED_THREAD_LEVEL MPI_THREAD_SINGLE
#endif

#ifndef TEST_MPI_RMA_ENABLE_CPE
#define TEST_MPI_RMA_ENABLE_CPE 0
#endif

#if TEST_MPI_RMA_ENABLE_CPE
#include <athread.h>
extern void SLAVE_FUN(test_mpi_rma_cpe_work)();
static __uncached struct cpe_progress_task cpe_progress_shared_task;
#endif

/*
 * Standalone MPI-3 RMA diagnostic.
 *
 * The default build uses only standard C/POSIX and MPI APIs. Defining
 * TEST_MPI_RMA_ENABLE_CPE adds an optional Sunway CPE progress experiment.
 */

enum {
    RESULT_FIELD_COUNT = 6
};

enum progress_mode {
    PROGRESS_TARGET_BUSY,
    PROGRESS_TARGET_IPROBE,
    PROGRESS_TARGET_BARRIER,
    PROGRESS_TARGET_CPE_IPROBE
};

struct options {
    double hold_seconds;
    double compute_seconds;
    int iterations;
    uint64_t cpe_iterations;
};

struct counter_window {
    unsigned long long *base;
    MPI_Win win;
    int lock_all_active;
};

struct progress_result {
    double fetch_seconds;
    double completion_seconds;
    double total_seconds;
};

struct stream_summary {
    double max_elapsed;
    double chunks_per_second;
};

static int world_rank;
static int world_size;
static volatile uint64_t busy_sink;

static const char *thread_level_name(int level)
{
    switch (level) {
    case MPI_THREAD_SINGLE:
        return "MPI_THREAD_SINGLE";
    case MPI_THREAD_FUNNELED:
        return "MPI_THREAD_FUNNELED";
    case MPI_THREAD_SERIALIZED:
        return "MPI_THREAD_SERIALIZED";
    case MPI_THREAD_MULTIPLE:
        return "MPI_THREAD_MULTIPLE";
    default:
        return "unknown";
    }
}

static void mpi_fail(int error_code, const char *expression,
                     const char *file, int line)
{
    char error_string[MPI_MAX_ERROR_STRING];
    int error_length = 0;

    MPI_Error_string(error_code, error_string, &error_length);
    fprintf(stderr,
            "[rank %d] MPI failure at %s:%d: %s: %.*s\n",
            world_rank, file, line, expression,
            error_length, error_string);
    MPI_Abort(MPI_COMM_WORLD, error_code);
    exit(EXIT_FAILURE);
}

#define MPI_CHECK(call)                                                      \
    do {                                                                     \
        int mpi_check_result_ = (call);                                      \
        if (mpi_check_result_ != MPI_SUCCESS)                                \
            mpi_fail(mpi_check_result_, #call, __FILE__, __LINE__);          \
    } while (0)

static double monotonic_seconds(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        perror("clock_gettime");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return (double)now.tv_sec + (double)now.tv_nsec * 1.0e-9;
}

/* Deliberately performs no MPI calls. */
static void busy_wait_no_mpi(double seconds)
{
    double deadline = monotonic_seconds() + seconds;
    uint64_t value = busy_sink + (uint64_t)(world_rank + 1);

    while (monotonic_seconds() < deadline) {
        value = value * UINT64_C(2862933555777941757) +
                UINT64_C(3037000493);
    }
    busy_sink = value;
}

static void target_progress_loop(double seconds)
{
    double deadline = monotonic_seconds() + seconds;

    while (monotonic_seconds() < deadline) {
        int flag;
        MPI_Status status;

        MPI_CHECK(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                             &flag, &status));
    }
}

#if TEST_MPI_RMA_ENABLE_CPE
static double target_progress_while_cpe_runs(uint64_t iterations,
                                             uint64_t *probe_count)
{
    double begin;

    memset((void *)&cpe_progress_shared_task, 0,
           sizeof(cpe_progress_shared_task));
    cpe_progress_shared_task.iterations = iterations;
    *probe_count = 0;
    begin = monotonic_seconds();
    __real_athread_spawn((void *)slave_test_mpi_rma_cpe_work,
                         (void *)&cpe_progress_shared_task, 1);

    while (cpe_progress_shared_task.done == 0) {
        int flag;
        MPI_Status status;

        MPI_CHECK(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                             &flag, &status));
        ++*probe_count;
    }
    athread_join();
    return monotonic_seconds() - begin;
}
#endif

static void initialize_local_counter(struct counter_window *window)
{
    unsigned long long zero = 0;
    unsigned long long previous = 0;

    MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, world_rank, 0, window->win));
    MPI_CHECK(MPI_Fetch_and_op(&zero, &previous, MPI_UNSIGNED_LONG_LONG,
                               world_rank, 0, MPI_REPLACE, window->win));
    MPI_CHECK(MPI_Win_unlock(world_rank, window->win));
}

static void create_distributed_counter_window(struct counter_window *window)
{
    memset(window, 0, sizeof(*window));
    MPI_CHECK(MPI_Win_allocate((MPI_Aint)sizeof(unsigned long long),
                               (int)sizeof(unsigned long long), MPI_INFO_NULL,
                               MPI_COMM_WORLD, &window->base, &window->win));
    initialize_local_counter(window);
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
}

static void create_central_counter_window(struct counter_window *window)
{
    MPI_Aint bytes = world_rank == 0
                   ? (MPI_Aint)sizeof(unsigned long long) : 0;

    memset(window, 0, sizeof(*window));
    MPI_CHECK(MPI_Win_allocate(bytes, (int)sizeof(unsigned long long),
                               MPI_INFO_NULL, MPI_COMM_WORLD,
                               &window->base, &window->win));
    MPI_CHECK(MPI_Win_lock_all(0, window->win));
    window->lock_all_active = 1;

    if (world_rank == 0) {
        unsigned long long zero = 0;
        unsigned long long previous = 0;

        MPI_CHECK(MPI_Fetch_and_op(&zero, &previous,
                                   MPI_UNSIGNED_LONG_LONG, 0, 0,
                                   MPI_REPLACE, window->win));
        MPI_CHECK(MPI_Win_flush(0, window->win));
    }
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
}

static void destroy_counter_window(struct counter_window *window)
{
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    if (window->lock_all_active)
        MPI_CHECK(MPI_Win_unlock_all(window->win));
    MPI_CHECK(MPI_Win_free(&window->win));
    memset(window, 0, sizeof(*window));
}

static const char *window_model_name(MPI_Win window)
{
    int *model = NULL;
    int flag = 0;

    MPI_CHECK(MPI_Win_get_attr(window, MPI_WIN_MODEL, &model, &flag));
    if (!flag || model == NULL) return "unknown";
    if (*model == MPI_WIN_UNIFIED) return "MPI_WIN_UNIFIED";
    if (*model == MPI_WIN_SEPARATE) return "MPI_WIN_SEPARATE";
    return "unknown";
}

static struct progress_result run_progress_case(
    const char *label, enum progress_mode mode,
    struct counter_window *output_window, double hold_seconds,
    uint64_t cpe_iterations)
{
    struct progress_result result = {0.0, 0.0, 0.0};
    double local_values[3] = {0.0, 0.0, 0.0};
    double root_values[3] = {0.0, 0.0, 0.0};
#if TEST_MPI_RMA_ENABLE_CPE
    double cpe_elapsed = 0.0;
    uint64_t cpe_probe_count = 0;
#else
    (void)cpe_iterations;
#endif
    double start_delay = hold_seconds < 0.5
                       ? hold_seconds * 0.1 : 0.05;

    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    if (world_rank == 0) {
        if (mode == PROGRESS_TARGET_BUSY) {
            busy_wait_no_mpi(hold_seconds);
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        } else if (mode == PROGRESS_TARGET_IPROBE) {
            target_progress_loop(hold_seconds);
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
#if TEST_MPI_RMA_ENABLE_CPE
        } else if (mode == PROGRESS_TARGET_CPE_IPROBE) {
            cpe_elapsed = target_progress_while_cpe_runs(
                cpe_iterations, &cpe_probe_count);
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
#endif
        } else {
            /* Stay inside MPI while rank 1 performs the RMA operation. */
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        }
    } else if (world_rank == 1) {
        unsigned long long increment = 1;
        unsigned long long previous = 0;
        double begin, after_fetch, after_completion;

        busy_wait_no_mpi(start_delay);
        begin = monotonic_seconds();
        MPI_CHECK(MPI_Fetch_and_op(&increment, &previous,
                                   MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_SUM,
                                   output_window->win));
        after_fetch = monotonic_seconds();
        MPI_CHECK(MPI_Win_flush(0, output_window->win));
        after_completion = monotonic_seconds();

        local_values[0] = after_fetch - begin;
        local_values[1] = after_completion - after_fetch;
        local_values[2] = after_completion - begin;
        MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    } else {
        if (mode == PROGRESS_TARGET_BARRIER) {
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        } else {
            busy_wait_no_mpi(hold_seconds);
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        }
    }

    MPI_CHECK(MPI_Reduce(local_values, root_values, 3, MPI_DOUBLE, MPI_MAX,
                         0, MPI_COMM_WORLD));
    if (world_rank == 0) {
        result.fetch_seconds = root_values[0];
        result.completion_seconds = root_values[1];
        result.total_seconds = root_values[2];
        printf("  %-34s fetch=%10.3f us  completion=%10.3f us"
               "  total=%10.3f us",
               label,
               result.fetch_seconds * 1.0e6,
               result.completion_seconds * 1.0e6,
               result.total_seconds * 1.0e6);
#if TEST_MPI_RMA_ENABLE_CPE
        if (mode == PROGRESS_TARGET_CPE_IPROBE)
            printf("  CPE=%8.3f s  probes=%" PRIu64,
                   cpe_elapsed, cpe_probe_count);
#endif
        putchar('\n');
    }
    return result;
}

static void print_progress_inference(struct progress_result busy,
                                     struct progress_result polling,
                                     struct progress_result barrier,
                                     double hold_seconds)
{
    double expected_wait = hold_seconds -
                         (hold_seconds < 0.5 ? hold_seconds * 0.1 : 0.05);
    double quick_threshold = expected_wait * 0.1;

    if (quick_threshold < 0.001) quick_threshold = 0.001;

    printf("\n  Interpretation\n");
    if (busy.completion_seconds > expected_wait * 0.5 &&
        (polling.completion_seconds < busy.completion_seconds * 0.25 ||
         barrier.completion_seconds < busy.completion_seconds * 0.25)) {
        printf("    RESULT: target-side MPI progress is required.\n"
               "    The RMA operation does not complete while rank 0 stays"
               " outside MPI,\n"
               "    but completes faster when rank 0 repeatedly enters MPI.\n"
               "    This strongly indicates a software/target-progress path"
               " rather than\n"
               "    fully asynchronous hardware offload for this operation.\n");
    } else if (busy.completion_seconds < quick_threshold) {
        printf("    RESULT: the operation progresses while rank 0 is outside"
               " MPI.\n"
               "    This is consistent with hardware offload or an MPI"
               " asynchronous\n"
               "    progress thread. Portable MPI cannot distinguish those"
               " two cases.\n");
    } else {
        printf("    RESULT: inconclusive. MPI_Iprobe may not drive RMA progress"
               " in this\n"
               "    implementation, or another transport bottleneck dominates.\n");
    }
}

static void input_ticket_operation(struct counter_window *input_window,
                                   int target, double *lock_seconds,
                                   double *fetch_seconds,
                                   double *unlock_seconds)
{
    unsigned long long increment = 1;
    unsigned long long ticket = 0;
    double begin, after_lock, after_fetch, after_unlock;

    begin = monotonic_seconds();
    MPI_CHECK(MPI_Win_lock(MPI_LOCK_SHARED, target, 0, input_window->win));
    after_lock = monotonic_seconds();
    MPI_CHECK(MPI_Fetch_and_op(&increment, &ticket,
                               MPI_UNSIGNED_LONG_LONG, target, 0, MPI_SUM,
                               input_window->win));
    after_fetch = monotonic_seconds();
    MPI_CHECK(MPI_Win_unlock(target, input_window->win));
    after_unlock = monotonic_seconds();

    *lock_seconds += after_lock - begin;
    *fetch_seconds += after_fetch - after_lock;
    *unlock_seconds += after_unlock - after_fetch;
}

static void output_offset_operation(struct counter_window *output_window,
                                    double *fetch_seconds,
                                    double *flush_seconds)
{
    const unsigned long long bytes = UINT64_C(64) << 20;
    unsigned long long offset = 0;
    double begin, after_fetch, after_flush;

    begin = monotonic_seconds();
    MPI_CHECK(MPI_Fetch_and_op(&bytes, &offset, MPI_UNSIGNED_LONG_LONG,
                               0, 0, MPI_SUM, output_window->win));
    after_fetch = monotonic_seconds();
    MPI_CHECK(MPI_Win_flush(0, output_window->win));
    after_flush = monotonic_seconds();

    *fetch_seconds += after_fetch - begin;
    *flush_seconds += after_flush - after_fetch;
}

static struct stream_summary run_stream_case(const char *label,
                                             int create_output_window,
                                             int use_output_rma,
                                             const struct options *options)
{
    struct counter_window input_window;
    struct counter_window output_window;
    struct stream_summary summary = {0.0, 0.0};
    double local[RESULT_FIELD_COUNT] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double *all = NULL;
    double begin;
    int iteration;

    create_distributed_counter_window(&input_window);
    memset(&output_window, 0, sizeof(output_window));
    if (create_output_window)
        create_central_counter_window(&output_window);

    if (world_rank == 0)
        all = (double *)calloc((size_t)world_size * RESULT_FIELD_COUNT,
                               sizeof(double));
    if (world_rank == 0 && all == NULL) {
        fprintf(stderr, "failed to allocate benchmark result buffer\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    begin = monotonic_seconds();
    for (iteration = 0; iteration < options->iterations; ++iteration) {
        int input_target = (world_rank + iteration + 1) % world_size;

        busy_wait_no_mpi(options->compute_seconds);
        input_ticket_operation(&input_window, input_target,
                               &local[1], &local[2], &local[3]);
        if (use_output_rma)
            output_offset_operation(&output_window, &local[4], &local[5]);
    }
    local[0] = monotonic_seconds() - begin;

    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    MPI_CHECK(MPI_Gather(local, RESULT_FIELD_COUNT, MPI_DOUBLE,
                         all, RESULT_FIELD_COUNT, MPI_DOUBLE,
                         0, MPI_COMM_WORLD));

    if (world_rank == 0) {
        double max_elapsed = 0.0;
        double input_sum = 0.0;
        double output_sum = 0.0;
        int rank;

        printf("\n  %s\n", label);
        printf("    rank    elapsed(s)   input-complete(us/op)"
               "   output-complete(us/op)\n");
        for (rank = 0; rank < world_size; ++rank) {
            const double *values = all + rank * RESULT_FIELD_COUNT;
            double input_total = values[1] + values[2] + values[3];
            double output_total = values[4] + values[5];

            if (values[0] > max_elapsed) max_elapsed = values[0];
            input_sum += input_total;
            output_sum += output_total;
            printf("    %4d    %10.6f       %14.3f"
                   "          %14.3f\n",
                   rank, values[0],
                   input_total * 1.0e6 / options->iterations,
                   use_output_rma
                       ? output_total * 1.0e6 / options->iterations : 0.0);
        }

        summary.max_elapsed = max_elapsed;
        summary.chunks_per_second =
            (double)world_size * options->iterations / max_elapsed;
        printf("    max rank elapsed             %12.6f s\n",
               max_elapsed);
        printf("    aggregate simulated chunks   %12d\n",
               world_size * options->iterations);
        printf("    aggregate chunk throughput   %12.3f chunks/s\n",
               summary.chunks_per_second);
        printf("    mean input RMA completion    %12.3f us/op\n",
               input_sum * 1.0e6 /
               ((double)world_size * options->iterations));
        if (use_output_rma) {
            printf("    mean output RMA completion   %12.3f us/op\n",
                   output_sum * 1.0e6 /
                   ((double)world_size * options->iterations));
        }
    }

    free(all);
    if (create_output_window)
        destroy_counter_window(&output_window);
    destroy_counter_window(&input_window);
    return summary;
}

static double parse_positive_double(const char *text, const char *option)
{
    char *end = NULL;
    double value;

    errno = 0;
    value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || value <= 0.0) {
        if (world_rank == 0)
            fprintf(stderr, "invalid value for %s: %s\n", option, text);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return value;
}

static int parse_positive_int(const char *text, const char *option)
{
    char *end = NULL;
    long value;

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        value <= 0 || value > INT_MAX) {
        if (world_rank == 0)
            fprintf(stderr, "invalid value for %s: %s\n", option, text);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static uint64_t parse_positive_uint64(const char *text, const char *option)
{
    char *end = NULL;
    unsigned long long value;

    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value == 0) {
        if (world_rank == 0)
            fprintf(stderr, "invalid value for %s: %s\n", option, text);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return (uint64_t)value;
}

static void parse_options(int argc, char **argv, struct options *options)
{
    int i;

    options->hold_seconds = 2.0;
    options->compute_seconds = 0.05;
    options->iterations = 20;
    options->cpe_iterations = UINT64_C(1000000000);

    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--hold") == 0 && i + 1 < argc) {
            options->hold_seconds =
                parse_positive_double(argv[++i], "--hold");
        } else if (strcmp(argv[i], "--compute-ms") == 0 && i + 1 < argc) {
            options->compute_seconds =
                parse_positive_double(argv[++i], "--compute-ms") / 1000.0;
        } else if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            options->iterations =
                parse_positive_int(argv[++i], "--iterations");
        } else if (strcmp(argv[i], "--cpe-iterations") == 0 &&
                   i + 1 < argc) {
            options->cpe_iterations =
                parse_positive_uint64(argv[++i], "--cpe-iterations");
        } else if (strcmp(argv[i], "--help") == 0) {
            if (world_rank == 0) {
                printf("Usage: %s [--hold SEC] [--compute-ms MS]"
                       " [--iterations N] [--cpe-iterations N]\n", argv[0]);
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        } else {
            if (world_rank == 0)
                fprintf(stderr, "unknown or incomplete option: %s\n", argv[i]);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }
    }
}

static void print_environment(MPI_Win probe_window,
                              const struct options *options,
                              int provided_thread_level)
{
    int mpi_major, mpi_minor;
    int library_length = 0;
    char library[MPI_MAX_LIBRARY_VERSION_STRING];

    MPI_CHECK(MPI_Get_version(&mpi_major, &mpi_minor));
    MPI_CHECK(MPI_Get_library_version(library, &library_length));

    printf("============================================================\n");
    printf("MPI RMA progress and dual-window diagnostic\n");
    printf("  ranks                 : %d\n", world_size);
    printf("  MPI standard          : %d.%d\n", mpi_major, mpi_minor);
    printf("  requested thread level: %s (%d)\n",
           thread_level_name(TEST_MPI_REQUIRED_THREAD_LEVEL),
           TEST_MPI_REQUIRED_THREAD_LEVEL);
    printf("  provided thread level : %s (%d)\n",
           thread_level_name(provided_thread_level),
           provided_thread_level);
    printf("  thread result         : %s\n",
           provided_thread_level >= TEST_MPI_REQUIRED_THREAD_LEVEL
               ? "sufficient" : "insufficient");
    printf("  RMA memory model      : %s\n", window_model_name(probe_window));
    printf("  target hold time      : %.3f s\n", options->hold_seconds);
    printf("  simulated compute     : %.3f ms/chunk\n",
           options->compute_seconds * 1000.0);
    printf("  iterations/rank       : %d\n", options->iterations);
#if TEST_MPI_RMA_ENABLE_CPE
    printf("  CPE progress test     : enabled\n");
    printf("  CPE iterations/core   : %" PRIu64 "\n",
           options->cpe_iterations);
#else
    printf("  CPE progress test     : disabled\n");
#endif
    printf("  MPI library           : %.*s\n",
           library_length, library);
    printf("============================================================\n");
}

int main(int argc, char **argv)
{
    struct options options;
    struct counter_window progress_window;
    struct progress_result busy_result;
    struct progress_result polling_result;
    struct progress_result barrier_result;
#if TEST_MPI_RMA_ENABLE_CPE
    struct progress_result cpe_polling_result;
#endif
    struct stream_summary input_only;
    struct stream_summary second_window_idle;
    struct stream_summary input_and_output;
    int provided_thread_level;

    if (sizeof(unsigned long long) != 8) {
        fprintf(stderr, "this test requires a 64-bit unsigned long long\n");
        return EXIT_FAILURE;
    }

    MPI_CHECK(MPI_Init_thread(&argc, &argv,
                              TEST_MPI_REQUIRED_THREAD_LEVEL,
                              &provided_thread_level));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
    MPI_CHECK(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN));

    if (world_size < 2) {
        if (world_rank == 0)
            fprintf(stderr, "run this test with at least two MPI ranks\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }

    parse_options(argc, argv, &options);
#if TEST_MPI_RMA_ENABLE_CPE
    if (world_rank == 0)
        athread_init();
#endif
    create_central_counter_window(&progress_window);
    if (world_rank == 0)
        print_environment(progress_window.win, &options,
                          provided_thread_level);

    if (world_rank == 0) {
        printf("\nTest 1: does rank 1 wait for rank 0 to enter MPI?\n");
        printf("  One Fetch_and_op targets rank 0; only completion time"
               " should change.\n");
    }
    busy_result = run_progress_case(
        "rank 0 busy, no MPI calls", PROGRESS_TARGET_BUSY,
        &progress_window, options.hold_seconds, options.cpe_iterations);
    polling_result = run_progress_case(
        "rank 0 polling MPI_Iprobe", PROGRESS_TARGET_IPROBE,
        &progress_window, options.hold_seconds, options.cpe_iterations);
    barrier_result = run_progress_case(
        "rank 0 waiting in MPI_Barrier", PROGRESS_TARGET_BARRIER,
        &progress_window, options.hold_seconds, options.cpe_iterations);
#if TEST_MPI_RMA_ENABLE_CPE
    cpe_polling_result = run_progress_case(
        "rank 0 CPE + MPI_Iprobe", PROGRESS_TARGET_CPE_IPROBE,
        &progress_window, options.hold_seconds, options.cpe_iterations);
#endif
    if (world_rank == 0)
        print_progress_inference(busy_result, polling_result,
                                 barrier_result, options.hold_seconds);
#if TEST_MPI_RMA_ENABLE_CPE
    if (world_rank == 0) {
        if (cpe_polling_result.completion_seconds <
            busy_result.completion_seconds * 0.25) {
            printf("    CPE overlap RESULT: polling MPI while the CPE runs"
                   " restores RMA progress.\n");
        } else {
            printf("    CPE overlap RESULT: MPI polling did not clearly"
                   " restore RMA progress.\n");
        }
    }
#endif
    destroy_counter_window(&progress_window);

    if (world_rank == 0) {
        printf("\nTest 2: input RMA versus input + output RMA\n");
        printf("  Input: distributed short lock/fetch/unlock ticket window.\n");
        printf("  Output: rank-0 lock_all/fetch/flush offset window.\n");
        printf("  No file I/O is performed; this isolates RMA behavior.\n");
    }
    input_only = run_stream_case(
        "A. input RMA only", 0, 0, &options);
    second_window_idle = run_stream_case(
        "B. input RMA + idle output window", 1, 0, &options);
    input_and_output = run_stream_case(
        "C. input RMA + active output RMA", 1, 1, &options);

    if (world_rank == 0) {
        printf("\n  Cross-case comparison\n");
        printf("    idle second-window slowdown   %12.3fx\n",
               second_window_idle.max_elapsed / input_only.max_elapsed);
        printf("    active input+output slowdown  %12.3fx\n",
               input_and_output.max_elapsed / input_only.max_elapsed);
        printf("    input-only throughput         %12.3f chunks/s\n",
               input_only.chunks_per_second);
        printf("    input+output throughput       %12.3f chunks/s\n",
               input_and_output.chunks_per_second);
        printf("\nNotes:\n");
        printf("  * MPI exposes no portable API that says whether an RMA"
               " operation used\n");
        printf("    a NIC atomic, shared-memory atomic, progress thread, or"
               " software message.\n");
        printf("  * Fast completion while rank 0 is outside MPI proves"
               " asynchronous\n");
        printf("    progress, but cannot distinguish hardware from a"
               " background thread.\n");
        printf("  * Completion tracking the hold time strongly indicates"
               " target-driven\n");
        printf("    software progress for this MPI/window/operation path.\n");
        printf("============================================================\n");
    }

    MPI_CHECK(MPI_Finalize());
    return EXIT_SUCCESS;
}

#endif /* TEST_MPI_RMA_BUILD_CPE_KERNEL */
