#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SWBWA_ALLOC_IMPLEMENTATION
#include "malloc_wrap.h"

int swbwa_test_cpe_id;
static int sdk_refuse;
static void *sdk_ptr;
static size_t sdk_bytes;
static uint32_t rng = 712398;

void *ldm_malloc(unsigned long bytes)
{
    void *p = sdk_refuse ? NULL : malloc(bytes);
    assert(!sdk_ptr);
    sdk_ptr = p; sdk_bytes = p ? bytes : 0;
    return p;
}
void ldm_free(void *p, unsigned long bytes)
{
    assert(p == sdk_ptr && bytes == sdk_bytes);
    free(p); sdk_ptr = NULL; sdk_bytes = 0;
}
void swbwa_cpe_fail(int code, long a, long b, long c)
{
    fprintf(stderr, "CPE FAIL %d %ld %ld %ld\n", code, a, b, c);
    abort();
}
static unsigned next(void)
{
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return rng;
}
static int is_ldm(void *p)
{
    return sdk_ptr && (uintptr_t)p >= (uintptr_t)sdk_ptr &&
           (uintptr_t)p - (uintptr_t)sdk_ptr < sdk_bytes;
}
static void check_bytes(const unsigned char *p, size_t size, unsigned char value)
{
    for (size_t i = 0; i < size; ++i) assert(p[i] == value);
}
extern int ksw_global2(int, const uint8_t *, int, const uint8_t *, int,
                      const int8_t *, int, int, int, int, int, int *, uint32_t **);
extern void *swbwa_test_unknown_malloc(int *bytes);
extern void *swbwa_test_unknown_calloc(int *bytes);
extern char *swbwa_test_unknown_strdup(const char *s);
extern void *swbwa_test_cached_query(void);
extern void swbwa_test_cached_query_free(void *ptr);

int main(void)
{
    size_t per_core = 32 << 20;
    char *pool = malloc(2 * per_core);
    uint64_t digest = UINT64_C(14695981039346656037);
    assert(pool);
    for (int core = 0; core < 2; ++core) {
        swbwa_test_cpe_id = core;
        set_big_buffer(pool, per_core);
        sdk_refuse = 0;
#if !SWBWA_CPE_MANUAL_LDM
        for (int site = 1; site <= 8; ++site)
            assert(swbwa_ldm_alloc(64, site) == NULL);
        assert(swbwa_ldm_outstanding() == 0 && swbwa_ldm_refusals() == 0);
#endif
        void *held = swbwa_ldm_alloc(SWBWA_LDM_SCRATCH_BUDGET_BYTES,
                                   SWBWA_LDM_AUTO_ARENA_SITE);
        assert(held);
        swbwa_ldm_allocator_begin();
        void *fallback = swbwa_auto_malloc(150, SWBWA_LDM_SITE_CHAIN);
        assert(!is_ldm(fallback));
        swbwa_auto_free(fallback);
        swbwa_ldm_allocator_end();
        assert(swbwa_ldm_outstanding() == SWBWA_LDM_SCRATCH_BUDGET_BYTES);
        swbwa_ldm_release(held, SWBWA_LDM_SCRATCH_BUDGET_BYTES);
        for (int refuse = 0; refuse < 2; ++refuse) {
            unsigned char *p[64] = {0};
            size_t sizes[64] = {0};
            sdk_refuse = refuse;
            swbwa_ldm_begin_batch();
            swbwa_ldm_allocator_begin();
#if SWBWA_CPE_LDM_ALLOC == 4
            void *cached_query = swbwa_test_cached_query();
            assert(!is_ldm(cached_query));
#endif
#if SWBWA_CPE_LDM_ALLOC >= 2 && SWBWA_CPE_LDM_BYTES > 4096
            if (!refuse) {
                /* Force spans across bitmap words, with live neighbours. */
                void *prefix = swbwa_auto_malloc(64, SWBWA_LDM_SITE_QUERY_DP);
                size_t half = SWBWA_CPE_LDM_BYTES / 2;
                unsigned char *wide = swbwa_auto_malloc(half, SWBWA_LDM_SITE_QUERY_DP);
                void *suffix = swbwa_auto_malloc(64, SWBWA_LDM_SITE_QUERY_DP);
                assert(is_ldm(prefix) && is_ldm(wide) && is_ldm(suffix));
                memset(prefix, 0x31, 64);
                memset(wide, 0x72, half);
                memset(suffix, 0x53, 64);
                wide = swbwa_auto_realloc(wide, half + 65, SWBWA_LDM_SITE_QUERY_DP);
                assert(!is_ldm(wide));
                check_bytes(wide, half, 0x72);
                check_bytes(prefix, 64, 0x31);
                check_bytes(suffix, 64, 0x53);
                swbwa_auto_free(prefix);
                swbwa_auto_free(wide);
                swbwa_auto_free(suffix);
            }
#endif
#if SWBWA_CPE_LDM_ALLOC >= 2 && SWBWA_CPE_LDM_BYTES > (128 << 10)
            if (!refuse) {
                /* Both an in-arena move and a heap spill must copy >64 KiB. */
                unsigned char *wide = swbwa_auto_malloc(70000, SWBWA_LDM_SITE_QUERY_DP);
                assert(is_ldm(wide));
                memset(wide, 0x69, 70000);
                wide = swbwa_auto_realloc(wide, 90000, SWBWA_LDM_SITE_QUERY_DP);
                assert(is_ldm(wide));
                check_bytes(wide, 70000, 0x69);
                memset(wide, 0x93, 90000);
                wide = swbwa_auto_realloc(wide, SWBWA_CPE_LDM_BYTES, SWBWA_LDM_SITE_QUERY_DP);
                assert(!is_ldm(wide));
                check_bytes(wide, 90000, 0x93);
                swbwa_auto_free(wide);
            }
#endif
            {
                int bytes = 75;
                void *q = swbwa_test_unknown_malloc(&bytes);
                assert(bytes == 76);
                assert(is_ldm(q) == (SWBWA_CPE_LDM_ALLOC == 3 && !refuse));
                swbwa_auto_free(q);
                q = swbwa_test_unknown_calloc(&bytes);
                assert(bytes == 77);
                assert(is_ldm(q) == (SWBWA_CPE_LDM_ALLOC == 3 && !refuse));
                check_bytes(q, 76, 0);
                swbwa_auto_free(q);
                swbwa_ldm_alloc_stats_t initial;
                swbwa_ldm_allocator_stats(&initial);
                assert(initial.site[SWBWA_LDM_SITE_OTHER].requests ==
                       (SWBWA_CPE_LDM_ALLOC == 2 || SWBWA_CPE_LDM_ALLOC == 4 ? 0 : 2));
#if SWBWA_CPE_LDM_ALLOC == 3
                q = swbwa_test_unknown_strdup("SAM across cross-runtime entries");
                assert(is_ldm(q) == !refuse);
                swbwa_ldm_allocator_suspend();
                swbwa_ldm_allocator_stats(&initial);
                assert(initial.carried_bytes == (refuse ? 0 : 64));
                swbwa_ldm_allocator_resume();
                assert(strcmp(q, "SAM across cross-runtime entries") == 0);
                swbwa_auto_free(q);
#endif
            }
#if SWBWA_CPE_LDM_ALLOC == 4
            assert(swbwa_ldm_alloc_site("slave_ksw_qinit_impl") == SWBWA_LDM_SITE_QUERY_DP);
            assert(swbwa_ldm_alloc_site("swbwa_extend2_scratch_reserve") == SWBWA_LDM_SITE_EXTEND_DP);
            assert(swbwa_ldm_alloc_site("slave_bwt_smem1a") == SWBWA_LDM_SITE_SMEM);
            assert(swbwa_ldm_alloc_site("smem_chain_alloc") == SWBWA_LDM_SITE_CHAIN_SEED);
            assert(swbwa_ldm_alloc_site("slave_mem_sort_dedup_patch") == SWBWA_LDM_SITE_DEDUP_SORT);
            assert(swbwa_ldm_alloc_tier(SWBWA_LDM_SITE_DEDUP_SORT) == 1);
            assert(swbwa_ldm_alloc_site("slave_mem_aln2sam") == SWBWA_LDM_SITE_OTHER);
            assert(swbwa_ldm_alloc_site("kputsn") == SWBWA_LDM_SITE_OTHER);
            if (!refuse && SWBWA_CPE_LDM_BYTES >= 4096) {
                size_t bytes = SWBWA_CPE_LDM_BYTES / 8;
                void *low[7];
                for (int i = 0; i < 7; ++i) {
                    low[i] = swbwa_auto_malloc(bytes, SWBWA_LDM_SITE_REFERENCE);
                    assert(is_ldm(low[i]) == (i < 6));
                }
                void *high = swbwa_auto_malloc(2 * bytes, SWBWA_LDM_SITE_QUERY_DP);
                assert(is_ldm(high));
                memset(high, 0x75, 2 * bytes);
                char *sam = swbwa_test_unknown_strdup("SAM stays in cross-segment heap");
                assert(!is_ldm(sam));
                for (int i = 0; i < 7; ++i) swbwa_auto_free(low[i]);
                check_bytes(high, 2 * bytes, 0x75);
                swbwa_auto_free(high);
                swbwa_ldm_allocator_suspend();
                swbwa_ldm_alloc_stats_t suspended;
                swbwa_ldm_allocator_stats(&suspended);
                assert(suspended.carried_bytes == 0 && suspended.reserved > 0);
                swbwa_ldm_allocator_resume();
                assert(strcmp(sam, "SAM stays in cross-segment heap") == 0);
                swbwa_auto_free(sam);
            }
#endif
            for (int i = 0; i < 16000; ++i) {
                unsigned slot = next() % 64;
                unsigned char value = (unsigned char)(slot + 17);
                size_t size = next() % 6000;
                check_bytes(p[slot], sizes[slot], value);
                if (!p[slot]) {
                    p[slot] = swbwa_auto_calloc(1, size, SWBWA_LDM_SITE_CHAIN);
                    check_bytes(p[slot], size, 0);
                } else if ((next() % 4) == 0) {
                    swbwa_auto_free(p[slot]); p[slot] = NULL; sizes[slot] = 0;
                    continue;
                } else {
                    p[slot] = swbwa_auto_realloc(p[slot], size, SWBWA_LDM_SITE_CHAIN);
                    check_bytes(p[slot], sizes[slot] < size ? sizes[slot] : size, value);
                }
                if (is_ldm(p[slot])) assert((uintptr_t)p[slot] % 64 == 0);
                if (size) memset(p[slot], value, size);
                sizes[slot] = size;
            }
            for (unsigned i = 0; i < 64; ++i) {
                check_bytes(p[i], sizes[i], (unsigned char)(i + 17));
                swbwa_auto_free(p[i]);
            }
            /* Beyond both LDM and pool size classes: preserve system fallback. */
            unsigned char *large = swbwa_auto_malloc(1 << 20, SWBWA_LDM_SITE_CHAIN);
            assert(!is_ldm(large));
            memset(large, 0x5a, 1 << 20);
            large = swbwa_auto_realloc(large, 2 << 20, SWBWA_LDM_SITE_CHAIN);
            assert(!is_ldm(large));
            check_bytes(large, 1 << 20, 0x5a);
            large = swbwa_auto_realloc(large, 75, SWBWA_LDM_SITE_CHAIN);
            assert(!is_ldm(large));
            check_bytes(large, 75, 0x5a);
            swbwa_auto_free(large);
            /* Only full mode admits unknown sites. Growth may spill to heap. */
            void *exported = swbwa_auto_malloc(75, SWBWA_LDM_SITE_OTHER);
            assert(is_ldm(exported) == (SWBWA_CPE_LDM_ALLOC == 3 && !refuse));
            exported = swbwa_auto_realloc(exported, 150, SWBWA_LDM_SITE_CHAIN);
            if (SWBWA_CPE_LDM_ALLOC != 3 || refuse) assert(!is_ldm(exported));
            swbwa_auto_free(exported);
            void *zero = swbwa_auto_malloc(0, SWBWA_LDM_SITE_CHAIN);
            swbwa_auto_free(zero);
            assert(swbwa_auto_realloc(NULL, 0, SWBWA_LDM_SITE_CHAIN) == NULL);
            swbwa_auto_free(NULL);

            for (int trial = 0; trial < 800; ++trial) {
                uint8_t q[300], t[360];
                int8_t mat[25];
                int qlen = 1 + next() % 300, tlen = 1 + next() % 360;
                int nc = 0, score;
                uint32_t *cigar = NULL;
                for (int i = 0; i < qlen; ++i) q[i] = next() % 5;
                for (int i = 0; i < tlen; ++i) t[i] = next() % 5;
                for (int i = 0; i < 25; ++i) mat[i] = i / 5 == i % 5 ? 2 : -4;
                /* A global traceback needs the terminal cell inside the band. */
                int band = abs(qlen - tlen) + 1 + next() % 50;
                score = ksw_global2(qlen, q, tlen, t, 5, mat, 6, 1, 6, 1, band, &nc, &cigar);
                printf("%d %d", score, nc);
                for (int i = 0; i < nc; ++i) printf(" %u", cigar[i]);
                putchar('\n');
                digest = (digest ^ (uint32_t)score) * UINT64_C(1099511628211);
                for (int i = 0; i < nc; ++i)
                    digest = (digest ^ cigar[i]) * UINT64_C(1099511628211);
                swbwa_auto_free(cigar);
            }
            swbwa_ldm_alloc_stats_t stats;
            swbwa_ldm_allocator_stats(&stats);
            assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].requests == 2400);
#if SWBWA_CPE_LDM_ALLOC >= 2
            if (!refuse) assert(stats.site[SWBWA_LDM_SITE_GLOBAL_DP].placed > 0);
#endif
            swbwa_ldm_allocator_end();
#if SWBWA_CPE_LDM_ALLOC == 4
            /* A public query survives destruction of the worker arena. */
            swbwa_test_cached_query_free(cached_query);
#endif
            assert(!sdk_ptr && swbwa_ldm_outstanding() == 0);
            assert(swbwa_ldm_peak() <= SWBWA_LDM_SCRATCH_BUDGET_BYTES);
        }
    }
    free(pool);
    printf("PASS allocator + 3200 exact global-DP cases: %016llx\n",
           (unsigned long long)digest);
    return 0;
}
