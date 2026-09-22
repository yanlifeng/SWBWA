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
} swbwa_ldm_site_stats_t;

typedef struct {
    swbwa_ldm_site_stats_t site[SWBWA_LDM_SITE_COUNT];
    unsigned long peak, misses, spills, arena_bytes, carried_bytes, reserved;
} swbwa_ldm_alloc_stats_t;

#if SWBWA_CPE_LDM_ALLOC
/* GCC folds the comparisons against __func__ at the allocation call site.
 * In auto mode these sites form the admission whitelist. In full mode they
 * are statistics labels only: unknown sites and SAM strings are eligible.
 * Tiered mode additionally labels audited worker-local objects. */
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

#ifdef __cplusplus
extern "C" {
#endif
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
