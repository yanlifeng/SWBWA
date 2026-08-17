#define _GNU_SOURCE

#include <stdint.h>

struct sleep_cpe_task {
    volatile uint64_t done;
    uint64_t iterations;
    volatile uint64_t checksum;
};

#ifdef TEST_USLEEP_BUILD_CPE

#include <slave.h>

void test_usleep_cpe_work(struct sleep_cpe_task *task)
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

#include <athread.h>
#include <mpi.h>

#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

extern void SLAVE_FUN(test_usleep_cpe_work)();

enum sleep_method {
    SLEEP_USLEEP_ONCE,
    SLEEP_NANOSLEEP_RETRY,
    SLEEP_CLOCK_NANOSLEEP_ABSOLUTE,
    SLEEP_METHOD_COUNT
};

struct sleep_stats {
    uint64_t logical_calls;
    uint64_t library_calls;
    uint64_t successful_calls;
    uint64_t interrupted_calls;
    uint64_t error_calls;
    uint64_t below_500us;
    uint64_t below_900us;
    double elapsed_seconds;
    double minimum_seconds;
    double maximum_seconds;
};

struct case_result {
    struct sleep_stats sleep;
    uint64_t probes;
    double wall_seconds;
};

static int world_rank;
static int world_size;
static __uncached struct sleep_cpe_task shared_cpe_task;

static void mpi_fail(int error_code, const char *expression,
                     const char *file, int line)
{
    char message[MPI_MAX_ERROR_STRING];
    int length = 0;

    MPI_Error_string(error_code, message, &length);
    fprintf(stderr, "[rank %d] MPI error at %s:%d: %s: %.*s\n",
            world_rank, file, line, expression, length, message);
    MPI_Abort(MPI_COMM_WORLD, error_code);
    exit(EXIT_FAILURE);
}

#define MPI_CHECK(call)                                                      \
    do {                                                                     \
        int mpi_result_ = (call);                                            \
        if (mpi_result_ != MPI_SUCCESS)                                      \
            mpi_fail(mpi_result_, #call, __FILE__, __LINE__);                \
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

static struct timespec monotonic_now(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        perror("clock_gettime");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return now;
}

static struct timespec add_nanoseconds(struct timespec value, long nanoseconds)
{
    value.tv_nsec += nanoseconds;
    while (value.tv_nsec >= 1000000000L) {
        value.tv_nsec -= 1000000000L;
        ++value.tv_sec;
    }
    return value;
}

static void initialize_sleep_stats(struct sleep_stats *stats)
{
    memset(stats, 0, sizeof(*stats));
    stats->minimum_seconds = DBL_MAX;
}

static const char *sleep_method_name(enum sleep_method method)
{
    switch (method) {
    case SLEEP_USLEEP_ONCE:
        return "usleep(1000), no retry";
    case SLEEP_NANOSLEEP_RETRY:
        return "nanosleep(1 ms), retry remaining time";
    case SLEEP_CLOCK_NANOSLEEP_ABSOLUTE:
        return "clock_nanosleep(absolute +1 ms), retry EINTR";
    default:
        return "unknown";
    }
}

static void record_sleep_result(struct sleep_stats *stats, double elapsed)
{
    ++stats->logical_calls;
    stats->elapsed_seconds += elapsed;
    if (elapsed < stats->minimum_seconds) stats->minimum_seconds = elapsed;
    if (elapsed > stats->maximum_seconds) stats->maximum_seconds = elapsed;
    if (elapsed < 0.0005) ++stats->below_500us;
    if (elapsed < 0.0009) ++stats->below_900us;
}

static void sleep_one_millisecond(enum sleep_method method,
                                  struct sleep_stats *stats)
{
    double begin = monotonic_seconds();
    int completed = 0;

    if (method == SLEEP_USLEEP_ONCE) {
        int result;
        int saved_errno;

        errno = 0;
        result = usleep(1000);
        saved_errno = errno;
        ++stats->library_calls;
        if (result == 0) {
            completed = 1;
        } else if (saved_errno == EINTR) {
            ++stats->interrupted_calls;
        } else {
            ++stats->error_calls;
        }
    } else if (method == SLEEP_NANOSLEEP_RETRY) {
        struct timespec request = {0, 1000000L};
        struct timespec remaining;

        for (;;) {
            int result;
            int saved_errno;

            errno = 0;
            result = nanosleep(&request, &remaining);
            saved_errno = errno;
            ++stats->library_calls;
            if (result == 0) {
                completed = 1;
                break;
            }
            if (saved_errno != EINTR) {
                ++stats->error_calls;
                break;
            }
            ++stats->interrupted_calls;
            request = remaining;
        }
    } else {
        struct timespec deadline = add_nanoseconds(monotonic_now(), 1000000L);

        for (;;) {
            int result = clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME,
                                         &deadline, NULL);

            ++stats->library_calls;
            if (result == 0) {
                completed = 1;
                break;
            }
            if (result != EINTR) {
                ++stats->error_calls;
                break;
            }
            ++stats->interrupted_calls;
        }
    }

    if (completed) ++stats->successful_calls;
    record_sleep_result(stats, monotonic_seconds() - begin);
}

static void mpi_probe(uint64_t *probe_count)
{
    int flag;
    MPI_Status status;

    MPI_CHECK(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                         &flag, &status));
    ++*probe_count;
}

static struct case_result run_fixed_case(enum sleep_method method,
                                         uint64_t iterations,
                                         int use_mpi_probe)
{
    struct case_result result;
    uint64_t i;
    double begin;

    memset(&result, 0, sizeof(result));
    initialize_sleep_stats(&result.sleep);
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    begin = monotonic_seconds();
    for (i = 0; i < iterations; ++i) {
        if (use_mpi_probe) mpi_probe(&result.probes);
        sleep_one_millisecond(method, &result.sleep);
    }
    result.wall_seconds = monotonic_seconds() - begin;
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    return result;
}

static struct case_result run_cpe_case(enum sleep_method method,
                                       uint64_t cpe_iterations,
                                       int use_mpi_probe)
{
    struct case_result result;
    double begin;

    memset(&result, 0, sizeof(result));
    initialize_sleep_stats(&result.sleep);
    memset((void *)&shared_cpe_task, 0, sizeof(shared_cpe_task));
    shared_cpe_task.iterations = cpe_iterations;

    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    begin = monotonic_seconds();
    __real_athread_spawn((void *)slave_test_usleep_cpe_work,
                         (void *)&shared_cpe_task, 1);
    while (shared_cpe_task.done == 0) {
        if (use_mpi_probe) mpi_probe(&result.probes);
        if (shared_cpe_task.done != 0) break;
        sleep_one_millisecond(method, &result.sleep);
    }
    athread_join();
    result.wall_seconds = monotonic_seconds() - begin;
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    return result;
}

static uint64_t minimum_u64(uint64_t current, uint64_t value)
{
    return value < current ? value : current;
}

static uint64_t maximum_u64(uint64_t current, uint64_t value)
{
    return value > current ? value : current;
}

static void print_case_summary(const char *context,
                               enum sleep_method method,
                               const struct case_result *local)
{
    struct case_result *all = NULL;
    uint64_t logical_min = UINT64_MAX;
    uint64_t logical_max = 0;
    uint64_t probe_min = UINT64_MAX;
    uint64_t probe_max = 0;
    uint64_t logical_total = 0;
    uint64_t library_total = 0;
    uint64_t success_total = 0;
    uint64_t interrupted_total = 0;
    uint64_t error_total = 0;
    uint64_t below_500us_total = 0;
    uint64_t below_900us_total = 0;
    double sleep_elapsed_total = 0.0;
    double sleep_minimum = DBL_MAX;
    double sleep_maximum = 0.0;
    double wall_minimum = DBL_MAX;
    double wall_maximum = 0.0;
    int rank;

    if (world_rank == 0) {
        all = calloc((size_t)world_size, sizeof(*all));
        if (all == NULL) {
            perror("calloc");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }
    }
    MPI_CHECK(MPI_Gather((void *)local, (int)sizeof(*local), MPI_BYTE,
                         all, (int)sizeof(*local), MPI_BYTE,
                         0, MPI_COMM_WORLD));
    if (world_rank != 0) return;

    for (rank = 0; rank < world_size; ++rank) {
        const struct case_result *value = &all[rank];

        logical_min = minimum_u64(logical_min, value->sleep.logical_calls);
        logical_max = maximum_u64(logical_max, value->sleep.logical_calls);
        probe_min = minimum_u64(probe_min, value->probes);
        probe_max = maximum_u64(probe_max, value->probes);
        logical_total += value->sleep.logical_calls;
        library_total += value->sleep.library_calls;
        success_total += value->sleep.successful_calls;
        interrupted_total += value->sleep.interrupted_calls;
        error_total += value->sleep.error_calls;
        below_500us_total += value->sleep.below_500us;
        below_900us_total += value->sleep.below_900us;
        sleep_elapsed_total += value->sleep.elapsed_seconds;
        if (value->sleep.logical_calls != 0 &&
            value->sleep.minimum_seconds < sleep_minimum)
            sleep_minimum = value->sleep.minimum_seconds;
        if (value->sleep.maximum_seconds > sleep_maximum)
            sleep_maximum = value->sleep.maximum_seconds;
        if (value->wall_seconds < wall_minimum)
            wall_minimum = value->wall_seconds;
        if (value->wall_seconds > wall_maximum)
            wall_maximum = value->wall_seconds;
    }

    printf("\n[%s]\n", context);
    printf("  method                     : %s\n", sleep_method_name(method));
    printf("  wall time across ranks     : %.6f - %.6f s\n",
           wall_minimum, wall_maximum);
    printf("  logical sleeps per rank    : %" PRIu64 " - %" PRIu64 "\n",
           logical_min, logical_max);
    printf("  MPI_Iprobe calls per rank  : %" PRIu64 " - %" PRIu64 "\n",
           probe_min, probe_max);
    printf("  logical/library calls      : %" PRIu64 " / %" PRIu64 "\n",
           logical_total, library_total);
    printf("  completed/EINTR/errors     : %" PRIu64 " / %" PRIu64
           " / %" PRIu64 "\n",
           success_total, interrupted_total, error_total);
    if (logical_total != 0) {
        printf("  measured sleep avg/min/max : %.3f / %.3f / %.3f us\n",
               sleep_elapsed_total * 1.0e6 / (double)logical_total,
               sleep_minimum * 1.0e6, sleep_maximum * 1.0e6);
        printf("  elapsed <500us / <900us    : %" PRIu64 " / %" PRIu64
               " (%.2f%% / %.2f%%)\n",
               below_500us_total, below_900us_total,
               100.0 * (double)below_500us_total / (double)logical_total,
               100.0 * (double)below_900us_total / (double)logical_total);
    }
    free(all);
}

static uint64_t parse_u64(const char *text, const char *name)
{
    char *end = NULL;
    unsigned long long value;

    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        if (world_rank == 0) fprintf(stderr, "invalid %s: %s\n", name, text);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }
    return (uint64_t)value;
}

int main(int argc, char **argv)
{
    uint64_t fixed_iterations = 500;
    uint64_t cpe_iterations = UINT64_C(500000000);
    int provided = MPI_THREAD_SINGLE;
    int method;

    MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
    if (provided < MPI_THREAD_MULTIPLE) {
        if (world_rank == 0)
            fprintf(stderr, "MPI_THREAD_MULTIPLE unavailable: provided=%d\n",
                    provided);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    if (argc > 1) fixed_iterations = parse_u64(argv[1], "fixed iterations");
    if (argc > 2) cpe_iterations = parse_u64(argv[2], "CPE iterations");

    athread_init();
    if (world_rank == 0) {
        printf("============================================================\n");
        printf("SWBWA usleep/MPI/CPE diagnostic\n");
        printf("  MPI ranks             : %d\n", world_size);
        printf("  MPI thread level      : %d (MULTIPLE=%d)\n",
               provided, MPI_THREAD_MULTIPLE);
        printf("  fixed sleep iterations: %" PRIu64 "\n", fixed_iterations);
        printf("  CPE loop iterations   : %" PRIu64 " per CPE\n", cpe_iterations);
        printf("============================================================\n");
    }

    for (method = 0; method < SLEEP_METHOD_COUNT; ++method) {
        struct case_result result = run_fixed_case(
            (enum sleep_method)method, fixed_iterations, 0);
        print_case_summary("fixed loop, sleep only",
                           (enum sleep_method)method, &result);
    }
    for (method = 0; method < SLEEP_METHOD_COUNT; ++method) {
        struct case_result result = run_fixed_case(
            (enum sleep_method)method, fixed_iterations, 1);
        print_case_summary("fixed loop, MPI_Iprobe then sleep",
                           (enum sleep_method)method, &result);
    }
    for (method = 0; method < SLEEP_METHOD_COUNT; ++method) {
        struct case_result result = run_cpe_case(
            (enum sleep_method)method, cpe_iterations, 0);
        print_case_summary("asynchronous CPE, sleep only",
                           (enum sleep_method)method, &result);
    }
    for (method = 0; method < SLEEP_METHOD_COUNT; ++method) {
        struct case_result result = run_cpe_case(
            (enum sleep_method)method, cpe_iterations, 1);
        print_case_summary("asynchronous CPE, MPI_Iprobe then sleep",
                           (enum sleep_method)method, &result);
    }

    MPI_CHECK(MPI_Finalize());
    return EXIT_SUCCESS;
}

#endif
