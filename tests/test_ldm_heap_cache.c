#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define SWBWA_ALLOC_IMPLEMENTATION
#include "malloc_wrap.h"

int swbwa_test_cpe_id;
static void *backing;
static size_t backing_bytes;
void *ldm_malloc(unsigned long bytes) {
    assert(!backing); backing = malloc(bytes); backing_bytes = bytes; return backing;
}
void ldm_free(void *p, unsigned long bytes) {
    assert(p == backing && bytes == backing_bytes); free(p); backing = NULL;
}
void swbwa_cpe_fail(int code, long a, long b, long c) {
    fprintf(stderr, "unexpected CPE error %d %ld %ld %ld\n", code, a, b, c); abort();
}

int main(void) {
    enum { N = 1200 };
    unsigned char *p[N];
    char *heap = malloc(128u << 20);
    assert(heap); set_big_buffer(heap, 128u << 20);
    swbwa_ldm_policy_t policy = {0};
    policy.scratch_hints = 255; policy.heap_cache = 1;
    policy.pool_bitmap = getenv("TEST_POOL_BITMAP") ? atoi(getenv("TEST_POOL_BITMAP")) : 0;
    swbwa_ldm_set_policy(&policy);
    for (int epoch = 0; epoch < 3; ++epoch) {
        swbwa_ldm_allocator_begin();
        /* More than 64 underlying segment trees: exercise a partial directory. */
        for (int i = 0; i < N; ++i) {
            p[i] = swbwa_auto_malloc(65536, SWBWA_LDM_SITE_OTHER);
            assert(p[i]); memset(p[i], i % 251, 256);
            p[i][65535] = i % 251;
        }
        /* Private payload stays valid while the temporary LDM directory dies. */
        swbwa_ldm_allocator_suspend(); swbwa_ldm_allocator_resume();
        for (int i = 0; i < N; i += 17) {
            p[i] = swbwa_auto_realloc(p[i], 65537, SWBWA_LDM_SITE_OTHER);
            assert(p[i][65535] == i % 251);
        }
        for (int i = 0; i < N; ++i) {
            int j = (i * 37) % N;
            for (int k = 0; k < 256; ++k) assert(p[j][k] == j % 251);
            assert(p[j][65535] == j % 251);
            swbwa_auto_free(p[j]);
        }
        swbwa_ldm_allocator_end();
        assert(!backing && swbwa_ldm_outstanding() == 0);
    }
    free(heap);
    puts("PASS heap lookup cache: >64 trees, gaps, realloc, reuse, suspend/resume");
    return 0;
}
