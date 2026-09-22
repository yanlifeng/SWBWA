#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <slave.h>
#include "swbwa_ldm_alloc.h"

#if SWBWA_CPE_LDM_ALLOC
#define SWBWA_ALLOC_IMPLEMENTATION
#include "malloc_wrap.h"

/* The original pool remains the sole fallback; its pointers/layout are unchanged. */
void *cpe_pool_malloc(size_t size);
void *cpe_pool_realloc(void *ptr, size_t size);
void cpe_pool_free(void *ptr);

#define BITMAP_WORDS ((SWBWA_CPE_LDM_BYTES / 64 + 63) / 64)
enum { UNIT = 64, SLOTS = SWBWA_CPE_LDM_BYTES / UNIT };
/* Preserve the small-arena layout; larger requests exceed a 16-bit byte count. */
#if SWBWA_CPE_LDM_BYTES > (96 << 10)
typedef uint32_t ldm_size_t;
#else
typedef uint16_t ldm_size_t;
#endif
typedef struct {
    void *backing;
    unsigned char *base;
    uint64_t used[BITMAP_WORDS];
    unsigned short spans[SLOTS];
    ldm_size_t sizes[SLOTS];
    unsigned live_bytes;
    unsigned suspended;
    swbwa_ldm_alloc_stats_t stats;
} ldm_allocator_t;

/* Indexed state follows the existing cross-runtime allocator convention.
 * Local handles (including opaque SAM pointers in the shared table) may only
 * be dereferenced/freed by the owning CPE. */
static ldm_allocator_t fallback_states[SWBWA_CPE_COUNT];
static ldm_allocator_t *allocators[SWBWA_CPE_COUNT];
static swbwa_ldm_alloc_stats_t completed_stats[SWBWA_CPE_COUNT];
#define BACKING_BYTES (sizeof(ldm_allocator_t) + SWBWA_CPE_LDM_BYTES + UNIT - 1)

static uint64_t slot_mask(unsigned count)
{
    return count == 64 ? UINT64_MAX : (UINT64_C(1) << count) - 1;
}

void swbwa_ldm_allocator_begin(void)
{
    ldm_allocator_t *a = &fallback_states[_MYID];
    if (allocators[_MYID])
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, 1, 0, 91);
#if SWBWA_CPE_LDM_ALLOC >= 2
    /* Hot bitmap/length metadata belongs in LDM too. Charge all metadata and
     * alignment slop against the SAME budget as manual LDM scratch. */
    void *backing = swbwa_ldm_alloc(BACKING_BYTES, SWBWA_LDM_AUTO_ARENA_SITE);
#if defined(__sw_slave__) && SWBWA_CPE_LDM_BYTES > (96 << 10)
    /* Large experimental arenas must leave room below the current LDM stack.
     * This is a refusal guard, not a proof of worst-case stack consumption. */
    volatile unsigned char stack_marker;
    uintptr_t stack = (uintptr_t)&stack_marker, begin = (uintptr_t)backing;
    unsigned long gap = backing && begin < stack && BACKING_BYTES <= stack - begin
                      ? stack - begin - BACKING_BYTES : 0;
    if (backing && gap < (32u << 10)) {
        swbwa_ldm_release(backing, BACKING_BYTES);
        backing = NULL;
    }
    if (_MYID == 0) {
        printf("[CPE LDM large] payload=%u backing=%lu begin=%#lx stack=%#lx gap=%lu accepted=%d\n",
               SWBWA_CPE_LDM_BYTES, (unsigned long)BACKING_BYTES,
               (unsigned long)begin, (unsigned long)stack, gap, backing != NULL);
    }
#endif
    if (backing != NULL) {
        a = backing;
        memset(a, 0, sizeof(*a));
        a->backing = backing;
        a->base = (unsigned char *)(((uintptr_t)(a + 1) + UNIT - 1) &
                                   ~(uintptr_t)(UNIT - 1));
        a->stats.arena_bytes = BACKING_BYTES;
    } else {
        memset(a, 0, sizeof(*a));
    }
#else
    memset(a, 0, sizeof(*a));
#endif
    allocators[_MYID] = a;
}

void swbwa_ldm_allocator_end(void)
{
    ldm_allocator_t *a = allocators[_MYID];
    /* Never hide an escaping allocation by resetting the arena. */
    if (!a || a->suspended)
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, 0, 0, 92);
    for (unsigned i = 0; i < BITMAP_WORDS; ++i)
        if (a->used[i])
            swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, i, (long)a->used[i], 92);
    completed_stats[_MYID] = a->stats;
    allocators[_MYID] = NULL;
    if (a->backing) swbwa_ldm_release(a->backing, BACKING_BYTES);
}

/* Keep the arena allocated while the host assigns final SAM destinations.
 * task_list preserves the CPE owner across the two cross-runtime entries. */
void swbwa_ldm_allocator_suspend(void)
{
    ldm_allocator_t *a = allocators[_MYID];
    if (!a || a->suspended)
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, 0, 0, 96);
#if SWBWA_CPE_LDM_ALLOC == 4
    if (a->live_bytes)
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, a->live_bytes, 0, 99);
#endif
    a->stats.carried_bytes = a->live_bytes;
    a->suspended = 1;
}

void swbwa_ldm_allocator_resume(void)
{
    ldm_allocator_t *a = allocators[_MYID];
    if (!a || !a->suspended)
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, 0, 0, 97);
    a->suspended = 0;
}

void swbwa_ldm_allocator_stats(swbwa_ldm_alloc_stats_t *stats)
{
    ldm_allocator_t *a = allocators[_MYID];
    *stats = a ? a->stats : completed_stats[_MYID];
}

static void mark_span(ldm_allocator_t *a, unsigned start, unsigned count, int occupied)
{
    while (count) {
        unsigned word = start / 64, bit = start % 64;
        unsigned take = count < 64 - bit ? count : 64 - bit;
        uint64_t mask = slot_mask(take) << bit;
        if (occupied) a->used[word] |= mask;
        else a->used[word] &= ~mask;
        start += take;
        count -= take;
    }
}

/* Return the first possible start past a conflicting bit, or the input start
 * when the complete span (which can cross several bitmap words) is free. */
static unsigned next_span(ldm_allocator_t *a, unsigned start, unsigned count)
{
    unsigned pos = start;
    while (count) {
        unsigned word = pos / 64, bit = pos % 64;
        unsigned take = count < 64 - bit ? count : 64 - bit;
        uint64_t conflict = a->used[word] & (slot_mask(take) << bit);
        if (conflict) return word * 64 + 64 - __builtin_clzll(conflict);
        pos += take;
        count -= take;
    }
    return start;
}

static void *try_alloc(ldm_allocator_t *a, size_t size, unsigned site)
{
    unsigned count, start = 0;
    if (!a) return NULL;
    if (a->suspended)
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, size, site, 98);
    if (site >= SWBWA_LDM_SITE_COUNT) site = SWBWA_LDM_SITE_OTHER;
    ++a->stats.site[site].requests;
    a->stats.site[site].bytes += size;
    /* Reserve room for several live objects; large requests never fragment it. */
#if SWBWA_CPE_LDM_ALLOC != 3
    if (site == SWBWA_LDM_SITE_OTHER) return NULL;
#endif
    if (size > SWBWA_CPE_LDM_BYTES / 2) return NULL;
    ++a->stats.site[site].small;
    if (!a->base) return NULL;
    count = (unsigned)((size ? size : 1) + UNIT - 1) / UNIT;
#if SWBWA_CPE_LDM_ALLOC == 4
    /* Lower-reuse buffers get at most 1/8 arena each and leave 1/4 free.
     * Global traceback is larger/less reused than EH/QP: apply the same
     * reserve to its >4 KiB requests, without changing the DP algorithm. */
    if (swbwa_ldm_alloc_tier(site) == 2 ||
        (site == SWBWA_LDM_SITE_GLOBAL_DP && size > 4096)) {
        if (count > SLOTS / 8 ||
            a->live_bytes + count * UNIT > (SLOTS - SLOTS / 4) * UNIT) {
            ++a->stats.reserved;
            return NULL;
        }
    }
#endif
    while (start + count <= SLOTS) {
        unsigned next = next_span(a, start, count);
        if (next == start) {
            mark_span(a, start, count, 1);
            a->spans[start] = count;
            a->sizes[start] = (ldm_size_t)size;
            a->live_bytes += count * UNIT;
            if (a->live_bytes > a->stats.peak) a->stats.peak = a->live_bytes;
            ++a->stats.site[site].placed;
            return a->base + start * UNIT;
        }
        start = next;
    }
    ++a->stats.misses;
    return NULL;
}

static int local_slot(ldm_allocator_t *a, const void *ptr)
{
    uintptr_t p = (uintptr_t)ptr, base;
    unsigned slot;
    if (!a) return -1;
    base = (uintptr_t)a->base;
    if (!a->base || p < base || p - base >= SWBWA_CPE_LDM_BYTES) return -1;
    slot = (unsigned)(p - base) / UNIT;
    if ((p - base) % UNIT || !a->spans[slot])
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, (long)ptr, slot, 93);
    return (int)slot;
}

static void release_slot(ldm_allocator_t *a, unsigned slot)
{
    unsigned count = a->spans[slot];
    mark_span(a, slot, count, 0);
    a->spans[slot] = 0;
    a->sizes[slot] = 0;
    a->live_bytes -= count * UNIT;
}

void *swbwa_auto_malloc(size_t size, unsigned site)
{
    void *ptr = try_alloc(allocators[_MYID], size, site);
    if (!ptr) ptr = cpe_pool_malloc(size);
    if (!ptr && size) swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)size, site, 94);
    return ptr;
}

void *swbwa_auto_calloc(size_t count, size_t size, unsigned site)
{
    void *ptr;
    if (size && count > SIZE_MAX / size)
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_SIZE_OVERFLOW, (long)count, (long)size, site);
    size *= count;
    ptr = swbwa_auto_malloc(size, site);
    if (size) memset(ptr, 0, size);
    return ptr;
}

void swbwa_auto_free(void *ptr)
{
    ldm_allocator_t *a = allocators[_MYID];
    int slot;
    if (!ptr) return;
    slot = local_slot(a, ptr);
    if (slot >= 0) release_slot(a, (unsigned)slot);
    else cpe_pool_free(ptr);
}

void *swbwa_auto_realloc(void *ptr, size_t size, unsigned site)
{
    ldm_allocator_t *a = allocators[_MYID];
    int slot;
    void *next;
    if (!ptr) return size ? swbwa_auto_malloc(size, site) : NULL;
    if (!size) { swbwa_auto_free(ptr); return NULL; }
    slot = local_slot(a, ptr);
    if (slot < 0) {
        /* Preserve existing heap ownership. New/realloc(NULL) objects are
         * admitted normally; an already spilled object stays on the heap. */
        next = cpe_pool_realloc(ptr, size);
    } else {
        if (size <= (size_t)a->spans[slot] * UNIT) {
            a->sizes[slot] = (ldm_size_t)size;
            return ptr;
        }
        next = try_alloc(a, size, site);
        if (!next) {
            next = cpe_pool_malloc(size);
            if (next) ++a->stats.spills;
        }
        if (next) {
            memcpy(next, ptr, a->sizes[slot]);
            release_slot(a, (unsigned)slot);
        }
    }
    if (!next) swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)size, site, 95);
    return next;
}

char *swbwa_auto_strdup(const char *s, unsigned site)
{
    size_t bytes = strlen(s) + 1;
    char *p = swbwa_auto_malloc(bytes, site);
    memcpy(p, s, bytes);
    return p;
}
#endif
