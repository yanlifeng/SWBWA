#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SWBWA_ALLOC_IMPLEMENTATION
#include "malloc_wrap.h"

int swbwa_test_cpe_id;
static int refuse;
static void *backing;
static size_t backing_bytes;
void *ldm_malloc(unsigned long n)
{
    assert(!backing);
    if (refuse) return NULL;
    backing = malloc(n);
    backing_bytes = n;
    return backing;
}
void ldm_free(void *p, unsigned long n)
{
    assert(p == backing && n == backing_bytes);
    free(p);
    backing = NULL;
}
void swbwa_cpe_fail(int code, long a, long b, long c)
{
    fprintf(stderr, "unexpected CPE error %d %ld %ld %ld\n", code, a, b, c);
    abort();
}
static void *heap(size_t n) { return malloc(n); }
static void drop(void *p) { free(p); }

/* Match the chain's actual key width and generated allocation-site names. */
typedef struct { int64_t pos; uint64_t rest[6]; } chain_key_t;
static int compare(chain_key_t a, chain_key_t b) { return (a.pos > b.pos) - (a.pos < b.pos); }
#define calloc(n,s) swbwa_auto_calloc(n,s,swbwa_ldm_alloc_site(__func__))
#define realloc(p,s) swbwa_auto_realloc(p,s,swbwa_ldm_alloc_site(__func__))
#define free(p) swbwa_auto_free(p)
#include "kbtree.h"
KBTREE_INIT(chn, chain_key_t, compare)

static void mem_chain(void)
{
    kbtree_t(chn) *tree = kb_init(chn, 512);
    assert(tree);
    for (int i = 0; i < 3001; ++i) {
        chain_key_t k = {0};
        k.pos = (i * 37) % 3001;
        k.rest[0] = (uint64_t)k.pos * 13;
        kb_putp(chn, tree, &k);
    }
    for (int i = 0; i < 3001; ++i) {
        chain_key_t k = {0}; k.pos = i;
        chain_key_t *found = kb_getp(chn, tree, &k);
        assert(found && found->rest[0] == (uint64_t)i * 13);
    }
    int visited = 0;
#define visit(p) do { assert((p)->pos == visited); ++visited; } while (0)
    __kb_traverse(chain_key_t, tree, visit);
#undef visit
    assert(visited == 3001);
    kb_destroy(chn, tree);
}

int main(void)
{
    char *pool = heap(32 << 20);
    assert(pool);
    set_big_buffer(pool, 32 << 20);
    for (refuse = 0; refuse < 2; ++refuse) {
        swbwa_ldm_policy_t p = {0};
        p.cap[SWBWA_LDM_SITE_CHAIN] = 4096;
        p.reserve[SWBWA_LDM_SITE_CHAIN] = 8192;
        p.profile = 1;
        p.fast_realloc = 1;
        p.scratch_hints = 255;
        p.heap_cache = 1;
        swbwa_ldm_set_policy(&p);
        swbwa_ldm_allocator_begin();
        for (int n = 0; n < 20; ++n) mem_chain();
        swbwa_ldm_alloc_stats_t s;
        swbwa_ldm_allocator_stats(&s);
        assert(s.site[SWBWA_LDM_SITE_CHAIN].requests > 0);
        assert((s.site[SWBWA_LDM_SITE_CHAIN].placed > 0) == !refuse);
        swbwa_ldm_allocator_suspend();
        swbwa_ldm_allocator_resume();
        swbwa_ldm_allocator_end();
        assert(!backing && swbwa_ldm_outstanding() == 0);
    }
    drop(pool);
    puts("PASS B-tree 120040 inserts/lookups with pool exhaustion, fallback and cleanup");
    return 0;
}
