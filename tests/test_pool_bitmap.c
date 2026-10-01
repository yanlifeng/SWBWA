#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#define SWBWA_ALLOC_IMPLEMENTATION
#include "../src/slave/malloc_wrap.c"

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
    fprintf(stderr, "CPE failure %d %ld %ld %ld\n", code, a, b, c); abort();
}

int main(void) {
    const int counts[] = {1, 8, 64, 128, 8192, 65536};
    char *heap = malloc(64u << 20);
    swbwa_ldm_policy_t policy = {0};
    policy.scratch_hints = 255;
    policy.pool_bitmap = 2;
    swbwa_ldm_set_policy(&policy);
    for (unsigned c = 0; c < sizeof(counts) / sizeof(counts[0]); ++c) {
        int n = counts[c];
        unsigned char *live = calloc(n, 1);
        set_big_buffer(heap, 64u << 20);
        swbwa_ldm_allocator_begin();
        bitmap_modes[0] = 0;
        struct swbwa_segment_tree *reference = build_segment_tree(16, n);
        bitmap_modes[0] = 2;
        struct swbwa_segment_tree *candidate = build_segment_tree(16, n);
        SWBWA_ALLOC_TREES[0][0] = reference;
        SWBWA_ALLOC_TREES[0][1] = candidate;
        SWBWA_ALLOC_TREE_COUNTS[0] = 2;
        SWBWA_ALLOC_INITIALIZED = 1;
        uint32_t rng = 17;
        for (int step = 0; step < 100000; ++step) {
            rng = rng * 1664525u + 1013904223u;
            int pos = (int)(rng % n) + 1;
            int old = segment_allocation_state(reference, pos);
            assert(old == segment_allocation_state(candidate, pos));
            assert(old == live[pos - 1]);
            int value = old ? -1 : 1;
            update_segment_tree(reference, pos, value);
            update_segment_tree(candidate, pos, value);
            live[pos - 1] = !old;
            assert(find_free_segment(reference) == find_free_segment(candidate));
            if (step % 7919 == 0) {
                swbwa_ldm_allocator_suspend();
                assert(candidate->bitmap == candidate->bitmap_home);
                swbwa_ldm_allocator_resume();
            }
        }
        for (int i = 1; i <= n; ++i) {
            if (live[i - 1]) continue;
            update_segment_tree(reference, i, 1);
            update_segment_tree(candidate, i, 1);
        }
        assert(find_free_segment(candidate) == -1);
        for (int i = n; i; --i) {
            update_segment_tree(reference, i, -1);
            update_segment_tree(candidate, i, -1);
            assert(find_free_segment(reference) == i);
            assert(find_free_segment(candidate) == i);
        }
        swbwa_ldm_allocator_end();
        assert(!backing);
        free(live);
    }
    free(heap);
    puts("PASS bitmap/tree differential: 600000 updates, lowest-free order, full/empty, LDM overflow, epoch flush");
    return 0;
}
