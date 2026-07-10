#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if SWBWA_ENABLE_HOST_MALLOC_WRAPPER
#  undef SWBWA_ENABLE_HOST_MALLOC_WRAPPER
#endif

#include "malloc_wrap.h"

typedef struct {
    unsigned long long malloc_calls;
    unsigned long long calloc_calls;
    unsigned long long realloc_calls;
    unsigned long long strdup_calls;
    unsigned long long free_calls;
    unsigned long long failed_calls;
    unsigned long long malloc_bytes;
    unsigned long long calloc_bytes;
    unsigned long long realloc_bytes;
    unsigned long long strdup_bytes;
    size_t largest_request;
    const char *largest_file;
    const char *largest_func;
    unsigned int largest_line;
} host_malloc_stats_t;

typedef struct {
    unsigned long long vm_peak_kib;
    unsigned long long vm_size_kib;
    unsigned long long vm_hwm_kib;
    unsigned long long vm_rss_kib;
    int available;
} process_memory_stats_t;

static host_malloc_stats_t host_malloc_stats;

static void update_largest_request(size_t bytes, const char *file,
                                   unsigned int line, const char *func)
{
    if (bytes <= host_malloc_stats.largest_request) return;

    host_malloc_stats.largest_request = bytes;
    host_malloc_stats.largest_file = file;
    host_malloc_stats.largest_line = line;
    host_malloc_stats.largest_func = func;
}

static void record_request(unsigned long long *calls,
                           unsigned long long *requested_bytes,
                           size_t bytes, const char *file,
                           unsigned int line, const char *func)
{
    ++*calls;
    *requested_bytes += (unsigned long long)bytes;
    update_largest_request(bytes, file, line, func);
}

static void record_failure(size_t bytes)
{
    if (bytes > 0) ++host_malloc_stats.failed_calls;
}

static void read_process_memory_stats(process_memory_stats_t *stats)
{
    char line[256];
    FILE *status;

    memset(stats, 0, sizeof(*stats));
    status = fopen("/proc/self/status", "r");
    if (status == NULL) return;

    while (fgets(line, sizeof(line), status) != NULL) {
        if (sscanf(line, "VmPeak: %llu kB", &stats->vm_peak_kib) == 1)
            continue;
        if (sscanf(line, "VmSize: %llu kB", &stats->vm_size_kib) == 1)
            continue;
        if (sscanf(line, "VmHWM: %llu kB", &stats->vm_hwm_kib) == 1)
            continue;
        if (sscanf(line, "VmRSS: %llu kB", &stats->vm_rss_kib) == 1)
            continue;
    }
    fclose(status);
    stats->available = 1;
}

static void print_bytes(const char *label, unsigned long long bytes)
{
    const double mib = (double)bytes / (1024.0 * 1024.0);
    const double gib = mib / 1024.0;

    fprintf(stderr, "    %-42s %10.3f GiB  (%12.3f MiB)\n",
            label, gib, mib);
}

void swbwa_host_malloc_stats_init(void)
{
    memset(&host_malloc_stats, 0, sizeof(host_malloc_stats));
}

void swbwa_host_malloc_stats_print(void)
{
    process_memory_stats_t memory;
    unsigned long long total_calls;
    unsigned long long total_bytes;

    read_process_memory_stats(&memory);
    total_calls = host_malloc_stats.malloc_calls +
                  host_malloc_stats.calloc_calls +
                  host_malloc_stats.realloc_calls +
                  host_malloc_stats.strdup_calls;
    total_bytes = host_malloc_stats.malloc_bytes +
                  host_malloc_stats.calloc_bytes +
                  host_malloc_stats.realloc_bytes +
                  host_malloc_stats.strdup_bytes;

    fprintf(stderr,
            "\n"
            "===================== Host Memory Report =====================\n"
            "  Process memory\n");
    if (memory.available) {
        print_bytes("peak virtual memory (VmPeak)", memory.vm_peak_kib << 10);
        print_bytes("current virtual memory (VmSize)", memory.vm_size_kib << 10);
        print_bytes("peak resident memory (VmHWM)", memory.vm_hwm_kib << 10);
        print_bytes("current resident memory (VmRSS)", memory.vm_rss_kib << 10);
    } else {
        fprintf(stderr, "    /proc/self/status is unavailable\n");
    }

    fprintf(stderr, "\n  Wrapped allocation requests\n");
    fprintf(stderr, "    %-42s %10llu\n", "malloc calls", host_malloc_stats.malloc_calls);
    fprintf(stderr, "    %-42s %10llu\n", "calloc calls", host_malloc_stats.calloc_calls);
    fprintf(stderr, "    %-42s %10llu\n", "realloc calls", host_malloc_stats.realloc_calls);
    fprintf(stderr, "    %-42s %10llu\n", "strdup calls", host_malloc_stats.strdup_calls);
    fprintf(stderr, "    %-42s %10llu\n", "free calls", host_malloc_stats.free_calls);
    fprintf(stderr, "    %-42s %10llu\n", "allocation calls", total_calls);
    fprintf(stderr, "    %-42s %10llu\n", "failed allocation calls", host_malloc_stats.failed_calls);

    fprintf(stderr, "\n  Requested bytes (cumulative, not concurrent)\n");
    print_bytes("malloc", host_malloc_stats.malloc_bytes);
    print_bytes("calloc", host_malloc_stats.calloc_bytes);
    print_bytes("realloc target sizes", host_malloc_stats.realloc_bytes);
    print_bytes("strdup", host_malloc_stats.strdup_bytes);
    print_bytes("total requested", total_bytes);
    print_bytes("largest single request", host_malloc_stats.largest_request);
    if (host_malloc_stats.largest_file != NULL) {
        fprintf(stderr, "    largest request site: %s:%u (%s)\n",
                host_malloc_stats.largest_file,
                host_malloc_stats.largest_line,
                host_malloc_stats.largest_func);
    }
    fprintf(stderr,
            "==============================================================\n"
            "\n");
}

void *wrap_malloc(size_t size, const char *file, unsigned int line, const char *func)
{
    void *ptr;

    record_request(&host_malloc_stats.malloc_calls,
                   &host_malloc_stats.malloc_bytes,
                   size, file, line, func);
    ptr = malloc(size);
    if (ptr == NULL) record_failure(size);
    return ptr;
}

void wrap_free(void *ptr, const char *file, unsigned int line, const char *func)
{
    (void)file;
    (void)line;
    (void)func;
    ++host_malloc_stats.free_calls;
    free(ptr);
}

void *wrap_realloc(void *ptr, size_t size, const char *file, unsigned int line, const char *func)
{
    void *new_ptr;

    record_request(&host_malloc_stats.realloc_calls,
                   &host_malloc_stats.realloc_bytes,
                   size, file, line, func);
    new_ptr = realloc(ptr, size);
    if (new_ptr == NULL && size > 0) record_failure(size);
    return new_ptr;
}

void *wrap_calloc(size_t nmemb, size_t size, const char *file, unsigned int line, const char *func)
{
    size_t bytes;
    void *ptr;

    if (size != 0 && nmemb > (size_t)-1 / size) {
        ++host_malloc_stats.calloc_calls;
        ++host_malloc_stats.failed_calls;
        return NULL;
    }

    bytes = nmemb * size;
    record_request(&host_malloc_stats.calloc_calls,
                   &host_malloc_stats.calloc_bytes,
                   bytes, file, line, func);
    ptr = malloc(bytes);
    if (ptr != NULL) {
        memset(ptr, 0, bytes);
    } else {
        record_failure(bytes);
    }
    return ptr;
}

char *wrap_strdup(const char *s, const char *file, unsigned int line, const char *func)
{
    size_t len;
    char *copy;

    if (s == NULL) {
        ++host_malloc_stats.strdup_calls;
        ++host_malloc_stats.failed_calls;
        return NULL;
    }

    len = strlen(s) + 1;
    record_request(&host_malloc_stats.strdup_calls,
                   &host_malloc_stats.strdup_bytes,
                   len, file, line, func);
    copy = malloc(len);
    if (copy != NULL) {
        memcpy(copy, s, len);
    } else {
        record_failure(len);
    }
    return copy;
}
