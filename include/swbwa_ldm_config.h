#ifndef SWBWA_LDM_CONFIG_H
#define SWBWA_LDM_CONFIG_H

/* Explicit private-object placement only: cache and stack are runtime settings. */
#define SWBWA_CPE_LDM_OFF    0
#define SWBWA_CPE_LDM_POOL   1
#define SWBWA_CPE_LDM_MANUAL 2

#ifndef SWBWA_CPE_LDM_MODE
#define SWBWA_CPE_LDM_MODE SWBWA_CPE_LDM_MANUAL
#endif
#if SWBWA_CPE_LDM_MODE < 0 || SWBWA_CPE_LDM_MODE > 2
#error "SWBWA_CPE_LDM_MODE must be 0 (off), 1 (pool), or 2 (manual)"
#endif

/* Derived implementation flags, not independent user options. Policy 4 is
 * the existing lifetime/access-tier allocator, including dedup scratch. */
#if defined(SWBWA_CPE_LDM_ALLOC) || defined(SWBWA_CPE_MANUAL_LDM)
#error "select SWBWA_CPE_LDM_MODE instead of legacy LDM switches"
#endif
#define SWBWA_CPE_MANUAL_LDM (SWBWA_CPE_LDM_MODE == SWBWA_CPE_LDM_MANUAL)
#define SWBWA_CPE_LDM_ALLOC (SWBWA_CPE_LDM_MODE == SWBWA_CPE_LDM_POOL ? 4 : 0)

#if SWBWA_CPE_LDM_MODE == SWBWA_CPE_LDM_POOL && \
    (SWBWA_USE_MPI || !SWBWA_USE_CROSS_SEGMENT || !SWBWA_ENABLE_CPE_MALLOC_WRAPPER)
#error "CPE LDM pool requires non-MPI cgs_cross+pool"
#endif

/* Internal capacity overrides are for tests/research, not another mode axis.
 * The ordinary pool is 32 KiB plus metadata, within the 40 KiB budget. */
#ifndef SWBWA_CPE_LDM_BYTES
#define SWBWA_CPE_LDM_BYTES (32 << 10)
#endif
#ifndef SWBWA_LDM_SCRATCH_BUDGET_BYTES
#define SWBWA_LDM_SCRATCH_BUDGET_BYTES (40 << 10)
#endif
#define SWBWA_LDM_AUTO_ARENA_SITE 9
#if SWBWA_CPE_LDM_BYTES < 256 || SWBWA_CPE_LDM_BYTES > (200 << 10) || SWBWA_CPE_LDM_BYTES % 64
#error "SWBWA_CPE_LDM_BYTES must be a multiple of 64 in [256, 204800]"
#endif
#if SWBWA_LDM_SCRATCH_BUDGET_BYTES < 0 || SWBWA_LDM_SCRATCH_BUDGET_BYTES > (224 << 10)
#error "tracked scratch budget must be in [0, 224 KiB]; reserve stack space separately"
#endif
/* Large research arenas need a smaller cache and per-batch stack validation. */
#if (SWBWA_CPE_LDM_BYTES > (96 << 10) || SWBWA_LDM_SCRATCH_BUDGET_BYTES > (128 << 10)) && SWBWA_CPE_MANUAL_LDM
#error "large LDM arena experiments require pool mode, without manual LDM"
#endif

#endif
