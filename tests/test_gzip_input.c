#include "swbwa_gzip.h"
#include "codec_job.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t check_bytes(FILE *expected, char *data, int64_t size)
{
    char buffer[8192];
    int64_t offset = 0;
    uint64_t lines = 0;
    assert(data && size > 0 && data[size] == 0 && data[size - 1] == '\n');
    while (offset < size) {
        size_t n = (size_t)(size - offset);
        if (n > sizeof(buffer)) n = sizeof(buffer);
        assert(fread(buffer, 1, n, expected) == n);
        assert(!memcmp(buffer, data + offset, n));
        offset += n;
    }
    for (offset = 0; offset < size; ++offset) lines += data[offset] == '\n';
    assert(lines % 4 == 0);
    return lines / 4;
}

int main(int argc, char **argv)
{
    swbwa_gzip_reader *reader;
    swbwa_gzip_metrics metrics;
    FILE *expected[2] = {NULL, NULL}, *report;
    char *held[2] = {NULL, NULL};
    int64_t held_size[2] = {0, 0};
    uint64_t total_records = 0, chunks = 0, bytes[2] = {0, 0};
    int expect_error, paired, failed = 0, layout;
    assert(argc == 12);
    paired = strcmp(argv[3], "-") != 0;
    expect_error = atoi(argv[6]);
    layout = atoi(argv[10]);
    expected[0] = fopen(argv[4], "rb");
    assert(expected[0]);
    if (paired) { expected[1] = fopen(argv[5], "rb"); assert(expected[1]); }
    swbwa_codec_init();
    reader = swbwa_gzip_open(argv[2], paired ? argv[3] : NULL,
        (size_t)strtoull(argv[7], NULL, 10), (unsigned)atoi(argv[8]), argv[1]);
    if (!reader) {
        assert(expect_error && *swbwa_gzip_error(NULL));
        fprintf(stderr, "expected open error: %s\n", swbwa_gzip_error(NULL));
        failed = 1;
    }
    if (reader) assert(!swbwa_gzip_set_batch_bytes(reader, (size_t)strtoull(argv[11], NULL, 10)));
    while (reader) {
        char *r1 = (char *)1, *r2 = (char *)1;
        int64_t n1 = -1, n2 = -1;
        int status = (chunks % 2 ? swbwa_gzip_next : swbwa_gzip_read_next)(reader, &r1, &n1, &r2, &n2);
        if (status < 0) {
            assert(expect_error && !r1 && !r2 && !n1 && !n2);
            fprintf(stderr, "expected read error: %s\n", swbwa_gzip_error(reader));
            failed = 1;
            assert(swbwa_gzip_read_next(reader, &r1, &n1, &r2, &n2) == -1);
            assert(!r1 && !r2 && !n1 && !n2);
            break;
        }
        if (!status) {
            assert(!r1 && !r2 && !n1 && !n2);
            assert(swbwa_gzip_read_next(reader, &r1, &n1, &r2, &n2) == 0);
            assert(!*swbwa_gzip_error(reader));
            break;
        }
        assert(status == 1);
        if (layout && swbwa_gzip_check_pe_layout(reader, r1, n1, r2, n2, (unsigned)layout)) {
            assert(expect_error);
            fprintf(stderr, "expected layout error: %s\n", swbwa_gzip_error(reader));
            free(r1); free(r2);
            failed = 1;
            assert(swbwa_gzip_read_next(reader, &r1, &n1, &r2, &n2) == -1);
            assert(!r1 && !r2 && !n1 && !n2);
            break;
        }
        {
            uint64_t records = check_bytes(expected[0], r1, n1);
            if (paired) assert(check_bytes(expected[1], r2, n2) == records);
            else assert(!r2 && !n2);
            total_records += records;
            ++chunks; bytes[0] += n1; bytes[1] += n2;
            assert(!swbwa_gzip_get_metrics(reader, &metrics));
            assert(metrics.records == total_records && metrics.last_records == records);
        }
        /* Model in-place formatting; retained allocations must survive next/close. */
        memset(r1, 0x5a, (size_t)n1);
        if (r2) memset(r2, 0x6b, (size_t)n2);
        if (!held[0]) {
            held[0] = r1; held[1] = r2; held_size[0] = n1; held_size[1] = n2;
        } else { free(r1); free(r2); }
    }
    assert(failed == expect_error);
    if (reader) {
        assert(!swbwa_gzip_get_metrics(reader, &metrics));
        if (!failed) {
            assert(metrics.chunks == chunks && metrics.records == total_records && !metrics.last_records);
            assert(total_records == strtoull(argv[9], NULL, 10));
            for (int i = 0; i <= paired; ++i) {
                assert(fgetc(expected[i]) == EOF);
                assert(metrics.input[i].decoded_bytes == bytes[i]);
                assert(metrics.input[i].delivered_bytes == bytes[i]);
                assert(metrics.input[i].copy_bytes == bytes[i]);
                assert(metrics.input[i].io_s >= 0 && metrics.input[i].prefetch_io_s >= 0);
                assert(metrics.input[i].discover_run_s >= 0 && metrics.input[i].resolve_run_s >= 0);
                assert(metrics.input[i].repair_s >= 0 && metrics.input[i].task_build_s >= 0);
                assert(metrics.input[i].input_prepare_s >= 0);
                assert(metrics.input[i].discover_wait_s >= 0 && metrics.input[i].resolve_wait_s >= 0);
                {
                    double runs = metrics.input[i].discover_run_s + metrics.input[i].resolve_run_s;
                    double waits = metrics.input[i].discover_wait_s + metrics.input[i].resolve_wait_s;
                    assert(runs >= metrics.input[i].kernel_s - 1e-8 &&
                           runs <= metrics.input[i].kernel_s + 1e-8);
                    assert(waits >= metrics.input[i].queue_s - 1e-8 &&
                           waits <= metrics.input[i].queue_s + 1e-8);
                }
            }
            printf("OK backend=%s paired=%d chunks=%llu records=%llu bytes=%llu jobs=%llu\n",
                argv[1], paired, (unsigned long long)chunks, (unsigned long long)total_records,
                (unsigned long long)bytes[0], (unsigned long long)metrics.input[0].jobs);
        }
        report = tmpfile(); assert(report);
        assert(!swbwa_gzip_report(reader, report));
        {
            char line[2048];
            rewind(report);
            assert(fgets(line, sizeof(line), report));
            assert(fgets(line, sizeof(line), report));
            assert(strstr(line, " io_s=") && strstr(line, " prefetch_io_s="));
            assert(strstr(line, " discover_run_s=") && strstr(line, " resolve_run_s="));
            assert(strstr(line, " repair_s=") && strstr(line, " input_prepare_s="));
            assert(strstr(line, " task_build_s=") && strstr(line, " discover_wait_s="));
            assert(strstr(line, " resolve_wait_s="));
        }
        fclose(report);
    }
    swbwa_gzip_close(reader);
    swbwa_gzip_close(NULL);
    swbwa_codec_finish();
    for (int i = 0; i < 2; ++i) {
        if (held[i]) {
            for (int64_t j = 0; j < held_size[i]; ++j) assert(held[i][j] == (i ? 0x6b : 0x5a));
            assert(!held[i][held_size[i]]);
            free(held[i]);
        }
        if (expected[i]) fclose(expected[i]);
    }
    return 0;
}
