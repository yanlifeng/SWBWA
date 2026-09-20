#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bwamem.h"
#include "swbwa_output.h"
#include "swbwa_discard_digest.h"
#include "swbwa_host_prep.h"

/* Only the members used by the extracted copy worker are needed here. */
typedef struct { const mem_opt_t *opt; bseq1_t *seqs; } worker_t;
static jmp_buf expected_failure;
static int failure_expected;

static void err_fatal(const char *function, const char *message, ...)
{
    if (failure_expected) longjmp(expected_failure, 1);
    fprintf(stderr, "unexpected failure in %s: %s\n", function, message);
    abort();
}

int swbwa_mpi_rank(void) { return 0; }
int swbwa_mpi_size(void) { return 1; }
int swbwa_mpi_is_root(void) { return 1; }
void swbwa_mpi_print_rank_ordered(void (*printer)(void)) { printer(); }

/* Independent scalar specification: one little-endian word per iteration. */
static uint64_t reference_hash(const void *data, size_t length)
{
    const unsigned char *p = data;
    const uint64_t m = UINT64_C(0xc6a4a7935bd1e995);
    uint64_t h = UINT64_C(0x9e3779b97f4a7c15) ^ ((uint64_t)length * m);
    size_t offset = 0;
    while (length - offset >= 8) {
        uint64_t word = 0;
        for (int i = 7; i >= 0; --i) word = (word << 8) | p[offset + (size_t)i];
        word *= m;
        word ^= word >> 47;
        word *= m;
        h ^= word;
        h *= m;
        offset += 8;
    }
    if (offset != length) {
        uint64_t tail = 0;
        for (size_t i = length; i > offset; --i) tail = (tail << 8) | p[i - 1];
        h ^= tail;
        h *= m;
    }
    h ^= h >> 47;
    h *= m;
    h ^= h >> 47;
    return h;
}

#include "actual_worker.inc"

static void allocate_slices(worker_t w, int n, int *lengths, char *block)
{
    struct { int *sam_lengths; } parameters = {lengths}, *para = &parameters;
    char *sam_blocks_static[SWBWA_PIPELINE_BUFFER_COUNT];
    unsigned long batch_number = 1;
    sam_blocks_static[0] = block;
#include "actual_slices.inc"
}

static void consume_records(bseq1_t *seqs, int n)
{
    struct { bseq1_t *seqs; int n_seqs; } batch = {seqs, n}, *data = &batch;
    int i;
#include "actual_output_loop.inc"
}

#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
static void check_hash_and_state(int enabled)
{
    unsigned char input[70016], metadata[32];
    for (size_t i = 0; i < sizeof(input); ++i) input[i] = (unsigned char)(i * 113 + 7);
    for (size_t offset = 0; offset < 8; ++offset) {
        for (size_t length = 0; length <= 1024; ++length)
            assert(swbwa_digest_hash_blob(input + offset, length) == reference_hash(input + offset, length));
        for (size_t length = 65529; length <= 65543; ++length)
            assert(swbwa_digest_hash_blob(input + offset, length) == reference_hash(input + offset, length));
    }
    for (size_t length = 1; length <= 65; ++length) {
        unsigned char *exact = malloc(length);
        assert(exact);
        memcpy(exact, input, length);
        assert(swbwa_digest_hash_blob(exact, length) == reference_hash(exact, length));
        free(exact);
    }
    uint64_t before = swbwa_digest_hash_blob(input, 70000);
    input[69999] ^= 1;
    assert(swbwa_digest_hash_blob(input, 70000) != before);
    assert(swbwa_digest_hash_blob(input, 70000) == reference_hash(input, 70000));
    assert(swbwa_output_write_digest(NULL) == -1 && errno == EINVAL);
    swbwa_digest_prepare(metadata + 1, 16, enabled);
    assert(swbwa_output_write_digest(metadata + 1) == -1);
    swbwa_digest_finish(metadata + 1, input, 15);
    assert(swbwa_output_write_digest(metadata + 1) == -1);
    swbwa_digest_prepare(metadata + 1, 16, !enabled);
    swbwa_digest_finish(metadata + 1, input, 16);
    assert(swbwa_output_write_digest(metadata + 1) == -1);
    swbwa_digest_prepare(metadata + 1, 16, enabled);
    swbwa_digest_finish(metadata + 1, NULL, 16);
    assert(swbwa_output_write_digest(metadata + 1) == -1);
    swbwa_digest_prepare(metadata + 1, 0, enabled);
    swbwa_digest_finish(metadata + 1, NULL, 0);
    assert(swbwa_output_write_digest(metadata + 1) == 0);
    if (!enabled) {
        swbwa_digest_prepare(metadata + 1, 70000, 0);
        swbwa_digest_finish(metadata + 1, (const void *)(uintptr_t)1, 70000);
        assert(swbwa_digest_load64(metadata + 17) == SWBWA_DIGEST_COUNT_READY);
    }
}
#endif

int main(int argc, char **argv)
{
    uint64_t expected_sum = 0, expected_xor = 0, expected_calls = 0, expected_bytes = 0;
    const char *policy = getenv("SWBWA_DISCARD_HASH");
    int enabled = !(policy && !strcmp(policy, "0"));
    assert(argc == 2);
    assert(swbwa_output_write(NULL, 0) == -1);
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
    assert(swbwa_output_discard_hash_active() == -1);
    assert(swbwa_output_write_digest(NULL) == -1);
#endif
    if (swbwa_output_open(argv[1], 0) != 0) return 9;
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
    assert(swbwa_output_discard_hash_active() == enabled);
    check_hash_and_state(enabled);
#endif
    for (int paired = 0; paired <= 1; ++paired) {
        mem_opt_t opt = {0};
        bseq1_t seqs[12];
        int lengths[12];
        char *blobs[12];
        worker_t w = { &opt, seqs };
        unsigned char block[SWBWA_CPE_FORMAT_BUFFER_BYTES + 16];
        opt.flag = paired ? MEM_F_PE : 0;
        for (int batch = 0; batch < 9; ++batch) {
            int n = batch == 0 ? 0 : (paired ? 12 : 11);
            memset(block, 0xa5, sizeof(block));
            memset(seqs, 0, sizeof(seqs));
            for (int i = 0; i < n; ++i) {
                lengths[i] = i == 0 ? 0 : (i * 31 + batch * 7) % 200;
                blobs[i] = malloc((size_t)lengths[i] + 1);
                assert(blobs[i]);
                for (int j = 0; j < lengths[i]; ++j)
                    blobs[i][j] = (j % 43 == 42) ? '\n' : (char)('!' + (i + j + batch) % 90);
                blobs[i][lengths[i]] = 0;
                if (lengths[i]) {
                    uint64_t hash = reference_hash(blobs[i], (size_t)lengths[i]);
                    ++expected_calls;
                    expected_bytes += (uint64_t)lengths[i];
                    if (enabled) { expected_sum += hash; expected_xor ^= hash; }
                }
            }
            allocate_slices(w, n, lengths, (char *)block);
            int tasks = paired ? n / 2 : n;
            int mid = tasks / 2;
            worker12_fast(&w, mid, tasks, 17, lengths, blobs, NULL);
            worker12_fast(&w, 0, mid, 3, lengths, blobs, NULL);
            for (int i = 0; i < n; ++i) {
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
                assert(seqs[i].sam == (char *)block + (size_t)i * SWBWA_DIGEST_RECORD_BYTES);
                assert(swbwa_digest_load64((unsigned char *)seqs[i].sam) == (uint64_t)lengths[i]);
#else
                assert(strlen(seqs[i].sam) == (size_t)lengths[i]);
#endif
            }
            consume_records(seqs, n);
            for (size_t i = SWBWA_CPE_FORMAT_BUFFER_BYTES; i < sizeof(block); ++i) assert(block[i] == 0xa5);
        }
        failure_expected = 1;
        lengths[0] = SWBWA_CPE_FORMAT_BUFFER_BYTES;
        if (setjmp(expected_failure) == 0) { allocate_slices(w, 1, lengths, (char *)block); abort(); }
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
        if (setjmp(expected_failure) == 0) {
            allocate_slices(w, SWBWA_CPE_FORMAT_BUFFER_BYTES / SWBWA_DIGEST_RECORD_BYTES + 1, lengths, (char *)block);
            abort();
        }
        lengths[0] = -1;
        if (setjmp(expected_failure) == 0) { allocate_slices(w, 1, lengths, (char *)block); abort(); }
#endif
        failure_expected = 0;
    }
    printf("calls=%" PRIu64 " bytes=%" PRIu64 " sum=0x%016" PRIx64 " xor=0x%016" PRIx64
           " enabled=%d hash_prefix_bytes=0\n", expected_calls, expected_bytes, expected_sum, expected_xor, enabled);
    assert(swbwa_output_flush() == 0);
    assert(swbwa_output_close() == 0);
    assert(swbwa_output_close() == 0);
#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
    assert(swbwa_output_write_digest(NULL) == -1);
#endif
    return 0;
}
