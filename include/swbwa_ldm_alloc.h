#ifndef SWBWA_LDM_ALLOC_H
#define SWBWA_LDM_ALLOC_H

#include <stddef.h>
#include <string.h>
#include "swbwa_config.h"

enum {
    SWBWA_LDM_SITE_OTHER,
    SWBWA_LDM_SITE_CHAIN,
    SWBWA_LDM_SITE_REFERENCE,
    SWBWA_LDM_SITE_GLOBAL_DP,
    SWBWA_LDM_SITE_REG2ALN,
    SWBWA_LDM_SITE_QUERY_DP,
    SWBWA_LDM_SITE_EXTEND_DP,
    SWBWA_LDM_SITE_SMEM,
    SWBWA_LDM_SITE_CHAIN_SEED,
    SWBWA_LDM_SITE_CONTEXT,
    SWBWA_LDM_SITE_DEDUP_SORT,
    SWBWA_LDM_SITE_COUNT
};

typedef struct {
    unsigned long requests, bytes, small, placed;
#if SWBWA_LDM_UNIFIED
    unsigned long histogram[6], placed_histogram[6];
    unsigned long peak_live, live, max_request, lifetime_events, releases;
#endif
} swbwa_ldm_site_stats_t;

typedef struct {
    swbwa_ldm_site_stats_t site[SWBWA_LDM_SITE_COUNT];
    unsigned long peak, misses, spills, arena_bytes, carried_bytes, reserved;
} swbwa_ldm_alloc_stats_t;

#if SWBWA_CPE_LDM_ALLOC
#ifdef __cplusplus
extern "C" {
#endif
#if SWBWA_LDM_UNIFIED
typedef struct {
    unsigned cap[SWBWA_LDM_SITE_COUNT];
    unsigned reserve[SWBWA_LDM_SITE_COUNT];
    unsigned profile;
    unsigned fast_realloc;
    unsigned scratch_hints;
    unsigned heap_cache;
    unsigned pool_bitmap;
} swbwa_ldm_policy_t;
/* Frozen production policy (B). Offline tuning is opt-in and never changes it. */
static inline void swbwa_ldm_policy_default(swbwa_ldm_policy_t *policy)
{
    const swbwa_ldm_policy_t defaults = {
        {0, 4096, 1024, 4096, 256, 0, 0, 0, 0, 0, 0},
        {0, 4096, 4096, 0, 4096, 0, 0, 0, 0, 0, 0},
        0, 1, 255, 1, 2
    };
    *policy = defaults;
}
void swbwa_pool_bitmap_begin(unsigned mode);
void swbwa_pool_bitmap_end(void);
void swbwa_ldm_set_policy(const swbwa_ldm_policy_t *policy);
void *swbwa_ldm_scratch_try(size_t size, unsigned legacy_site);
int swbwa_ldm_scratch_release(void *ptr);
#endif
/* GCC folds comparisons against __func__ at the allocation call site.
 * Only audited owner-local objects are eligible; unknown sites and SAM
 * strings stay on the heap. Scratch hints use a separate explicit path. */
static inline __attribute__((always_inline)) unsigned swbwa_ldm_alloc_site(const char *func)
{
    if (__builtin_strcmp(func, "mem_chain2aln") == 0 ||
        __builtin_strcmp(func, "slave_mem_chain2aln") == 0) return SWBWA_LDM_SITE_CHAIN;
    if (__builtin_strcmp(func, "bns_get_seq") == 0 ||
        __builtin_strcmp(func, "slave_bns_get_seq") == 0) return SWBWA_LDM_SITE_REFERENCE;
    if (__builtin_strcmp(func, "ksw_global2") == 0 ||
        __builtin_strcmp(func, "slave_ksw_global2") == 0) return SWBWA_LDM_SITE_GLOBAL_DP;
    if (__builtin_strcmp(func, "mem_reg2aln") == 0 ||
        __builtin_strcmp(func, "slave_mem_reg2aln") == 0) return SWBWA_LDM_SITE_REG2ALN;
#if SWBWA_CPE_LDM_ALLOC == 4
#define SWBWA_LDM_FUNC(name) (__builtin_strcmp(func, #name) == 0 || \
                            __builtin_strcmp(func, "slave_" #name) == 0)
    if (SWBWA_LDM_FUNC(ksw_qinit_impl) || SWBWA_LDM_FUNC(ksw_qinit_u8_pair))
        return SWBWA_LDM_SITE_QUERY_DP;
    if (SWBWA_LDM_FUNC(swbwa_extend2_scratch_reserve))
        return SWBWA_LDM_SITE_EXTEND_DP;
    if (SWBWA_LDM_FUNC(bwt_smem1a) || SWBWA_LDM_FUNC(mem_collect_intv))
        return SWBWA_LDM_SITE_SMEM;
    if (SWBWA_LDM_FUNC(smem_chain_alloc) || SWBWA_LDM_FUNC(test_and_merge))
        return SWBWA_LDM_SITE_CHAIN_SEED;
    if (SWBWA_LDM_FUNC(worker12_context_init) || SWBWA_LDM_FUNC(smem_aux_init))
        return SWBWA_LDM_SITE_CONTEXT;
    if (SWBWA_LDM_FUNC(mem_sort_dedup_patch))
        return SWBWA_LDM_SITE_DEDUP_SORT;
#if SWBWA_LDM_UNIFIED
    /* These containers are consumed and destroyed by the same CPE during
     * the first worker pass. The chain result stores values, not B-tree
     * node pointers. SAM strings and returned CIGAR storage stay excluded. */
    if (SWBWA_LDM_FUNC(kb_init_chn) || SWBWA_LDM_FUNC(__kb_split_chn) ||
        SWBWA_LDM_FUNC(kb_putp_chn) || SWBWA_LDM_FUNC(kb_destroy_chn) ||
        SWBWA_LDM_FUNC(mem_chain) || SWBWA_LDM_FUNC(mem_chain_flt) ||
        SWBWA_LDM_FUNC(mem_pair) || SWBWA_LDM_FUNC(mem_mark_primary_se_core) ||
        SWBWA_LDM_FUNC(mem_mark_primary_se) || SWBWA_LDM_FUNC(mem_sam_pe) ||
        SWBWA_LDM_FUNC(mem_reg2sam))
        return SWBWA_LDM_SITE_CHAIN;
    if (SWBWA_LDM_FUNC(swbwa_matesw_prepare))
        return SWBWA_LDM_SITE_REFERENCE;
#endif
#undef SWBWA_LDM_FUNC
#endif
    return SWBWA_LDM_SITE_OTHER;
}

/* Tier 0 stays on the heap. Tier 1 is recurrence-heavy DP/small metadata;
 * tier 2 is useful scratch, but must leave space for tier 1. */
static inline unsigned swbwa_ldm_alloc_tier(unsigned site)
{
    switch (site) {
    case SWBWA_LDM_SITE_GLOBAL_DP:
    case SWBWA_LDM_SITE_QUERY_DP:
    case SWBWA_LDM_SITE_EXTEND_DP:
    case SWBWA_LDM_SITE_CONTEXT:
    case SWBWA_LDM_SITE_DEDUP_SORT:
        return 1;
    case SWBWA_LDM_SITE_CHAIN:
    case SWBWA_LDM_SITE_REFERENCE:
    case SWBWA_LDM_SITE_REG2ALN:
    case SWBWA_LDM_SITE_SMEM:
    case SWBWA_LDM_SITE_CHAIN_SEED:
        return 2;
    default:
        return 0;
    }
}

void swbwa_ldm_allocator_begin(void);
void swbwa_ldm_allocator_end(void);
void swbwa_ldm_allocator_suspend(void);
void swbwa_ldm_allocator_resume(void);
void swbwa_ldm_allocator_stats(swbwa_ldm_alloc_stats_t *stats);
void *swbwa_auto_malloc(size_t size, unsigned site);
void *swbwa_auto_calloc(size_t count, size_t size, unsigned site);
void *swbwa_auto_realloc(void *ptr, size_t size, unsigned site);
char *swbwa_auto_strdup(const char *s, unsigned site);
void swbwa_auto_free(void *ptr);
#ifdef __cplusplus
}
#endif
#endif
#endif
