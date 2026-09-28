#include "swbwa_input.h"
#include "swbwa_host_workers.h"

#include <errno.h>
#include <limits.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

int swbwa_input_reader_count(const char *value)
{
    if (value == NULL)
        return SWBWA_HOST_MPE_THREADS > 1 ? SWBWA_HOST_MPE_THREADS : 0;
    if (value[0] < '0' || value[0] > '6' || value[1] != '\0' ||
        value[0] - '0' > SWBWA_HOST_MPE_THREADS) {
        errno = EINVAL;
        return -1;
    }
    return value[0] - '0';
}

static double input_now(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
}

void swbwa_input_usage_sample(swbwa_input_usage_t *usage)
{
    memset(usage, 0, sizeof(*usage));
#ifdef RUSAGE_THREAD
    {
        struct rusage ru;
        if (getrusage(RUSAGE_THREAD, &ru) != 0) return;
        usage->samples = 1;
        usage->user_seconds = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6;
        usage->system_seconds = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
        usage->minor_faults = ru.ru_minflt;
        usage->major_faults = ru.ru_majflt;
        usage->voluntary_switches = ru.ru_nvcsw;
        usage->involuntary_switches = ru.ru_nivcsw;
    }
#endif
}

void swbwa_input_usage_add(swbwa_input_usage_t *total,
                          const swbwa_input_usage_t *before,
                          const swbwa_input_usage_t *after)
{
    if (!before->samples || !after->samples) return;
    ++total->samples;
    total->user_seconds += after->user_seconds - before->user_seconds;
    total->system_seconds += after->system_seconds - before->system_seconds;
    total->minor_faults += after->minor_faults - before->minor_faults;
    total->major_faults += after->major_faults - before->major_faults;
    total->voluntary_switches += after->voluntary_switches - before->voluntary_switches;
    total->involuntary_switches += after->involuntary_switches - before->involuntary_switches;
}

typedef struct {
    const swbwa_input_slice_t *inputs;
    int input_count, readers, profile;
    swbwa_input_result_t results[SWBWA_HOST_MPE_THREADS];
    int errors[SWBWA_HOST_MPE_THREADS];
} input_task_t;

static size_t slice_boundary(size_t length, int index, int count)
{
    size_t remainder = length % (size_t)count;
    return length / (size_t)count * (size_t)index +
           ((size_t)index < remainder ? (size_t)index : remainder);
}

static void read_slices(void *opaque, int worker, int workers)
{
    input_task_t *task = opaque;
    swbwa_input_result_t *result;
    swbwa_input_usage_t before = {0}, after;
    int input;

    if (workers < task->readers) {
        task->errors[worker] = EINVAL;
        return;
    }
    if (worker >= task->readers) return;
    result = &task->results[worker];
    if (task->profile) swbwa_input_usage_sample(&before);
    for (input = 0; input < task->input_count; ++input) {
        const swbwa_input_slice_t *slice = &task->inputs[input];
        size_t pos = slice_boundary(slice->length, worker, task->readers);
        size_t end = slice_boundary(slice->length, worker + 1, task->readers);
        while (pos < end) {
            size_t length = end - pos > INT_MAX ? INT_MAX : end - pos;
            int64_t offset = slice->offset + (int64_t)pos;
            double start = input_now();
            ssize_t count = pread(slice->fd, slice->buffer + pos, length, (off_t)offset);
            int error = errno;
            double seconds = input_now() - start;
            ++result->calls;
            result->syscall_seconds += seconds;
            if (seconds > result->max_seconds) {
                result->max_seconds = seconds;
                result->max_input = input + 1;
                result->max_offset = offset;
                result->max_bytes = count > 0 ? (uint64_t)count : 0;
            }
            if (count < 0 && error == EINTR) continue;
            if (count <= 0) {
                task->errors[worker] = count == 0 ? EIO : error;
                goto done;
            }
            pos += (size_t)count;
            result->bytes += (uint64_t)count;
        }
    }
done:
    if (task->profile) {
        swbwa_input_usage_sample(&after);
        swbwa_input_usage_add(&result->usage, &before, &after);
    }
}

int swbwa_input_read(const swbwa_input_slice_t *inputs, int input_count,
                     int readers, int profile, swbwa_input_result_t *result)
{
    input_task_t task;
    double start;
    int i, error = 0;

    if (!inputs || !result || input_count < 1 || input_count > 2 ||
        readers < 1 || readers > SWBWA_HOST_MPE_THREADS) {
        errno = EINVAL;
        return -1;
    }
    for (i = 0; i < input_count; ++i) {
        if (inputs[i].offset < 0 ||
            (inputs[i].length && (!inputs[i].buffer || inputs[i].fd < 0))) {
            errno = EINVAL;
            return -1;
        }
        if ((uint64_t)inputs[i].length > (uint64_t)INT64_MAX - inputs[i].offset) {
            errno = EOVERFLOW;
            return -1;
        }
    }
    memset(result, 0, sizeof(*result));
    memset(&task, 0, sizeof(task));
    task.inputs = inputs;
    task.input_count = input_count;
    task.readers = readers;
    task.profile = profile;
    start = input_now();
    if (readers == 1) read_slices(&task, 0, 1);
    else swbwa_host_workers_run(read_slices, &task);
    result->wall_seconds = input_now() - start;
    for (i = 0; i < readers; ++i) {
        const swbwa_input_result_t *part = &task.results[i];
        const swbwa_input_usage_t *usage = &part->usage;
        if (task.errors[i] && !error) error = task.errors[i];
        result->bytes += part->bytes;
        result->calls += part->calls;
        result->syscall_seconds += part->syscall_seconds;
        if (part->max_seconds > result->max_seconds) {
            result->max_seconds = part->max_seconds;
            result->max_input = part->max_input;
            result->max_offset = part->max_offset;
            result->max_bytes = part->max_bytes;
        }
        result->usage.samples += usage->samples;
        result->usage.user_seconds += usage->user_seconds;
        result->usage.system_seconds += usage->system_seconds;
        result->usage.minor_faults += usage->minor_faults;
        result->usage.major_faults += usage->major_faults;
        result->usage.voluntary_switches += usage->voluntary_switches;
        result->usage.involuntary_switches += usage->involuntary_switches;
    }
    if (error) {
        errno = error;
        return -1;
    }
    return 0;
}
