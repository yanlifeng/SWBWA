#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include <slave.h>

#if SWBWA_ENABLE_CPE_MALLOC_WRAPPER
/* Don't wrap ourselves */
#  undef SWBWA_ENABLE_CPE_MALLOC_WRAPPER
#endif

#include "malloc_wrap.h"

#if SWBWA_LDM_UNIFIED
struct pool_bitmap {
    unsigned used, first_word;
    uint64_t words[];
};
#endif

struct swbwa_segment_tree {
    int *tree;
    int seg_size;
    int leaf_offset;
    char *start_address;
    char *end_address;
#if SWBWA_LDM_UNIFIED
    struct pool_bitmap *bitmap;
    struct pool_bitmap *bitmap_home;
#endif
};

enum {
    SWBWA_ALLOC_SIZE_CLASS_COUNT = 17,
    SWBWA_ALLOC_MAX_TREES_PER_CLASS = 1000
};

static const int block_sizes[SWBWA_ALLOC_SIZE_CLASS_COUNT] = {
    1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
    1 << 8, 1 << 9, 1 << 10, 1 << 11, 1 << 12, 1 << 13,
    1 << 14, 1 << 15, 1 << 16, 1 << 17, 1 << 18
};

static const int initial_tree_sizes[SWBWA_ALLOC_SIZE_CLASS_COUNT] = {
    8, 8, 8, 32, 8, 1024, 32, 8192, 128, 8, 8, 8, 8, 8, 8, 8, 8
};

int allocator_initialized[SWBWA_CPE_COUNT] = {0};
struct swbwa_segment_tree
    *allocator_trees[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT]
                    [SWBWA_ALLOC_MAX_TREES_PER_CLASS];
int tree_counts[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT];
int next_tree_sizes[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT];

#define SWBWA_ALLOC_INITIALIZED allocator_initialized[_MYID]
#define SWBWA_ALLOC_TREES allocator_trees[_MYID]
#define SWBWA_ALLOC_TREE_COUNTS tree_counts[_MYID]
#define SWBWA_ALLOC_NEXT_TREE_SIZES next_tree_sizes[_MYID]

static char *pool_starts[SWBWA_CPE_COUNT];
static size_t pool_offsets[SWBWA_CPE_COUNT];
static size_t pool_sizes[SWBWA_CPE_COUNT];

#if SWBWA_LDM_UNIFIED
static unsigned bitmap_modes[SWBWA_CPE_COUNT];
enum { POOL_BITMAP_CACHE_BYTES = 4096 };
static unsigned char *bitmap_caches[SWBWA_CPE_COUNT];
static unsigned bitmap_offsets[SWBWA_CPE_COUNT];

static size_t bitmap_bytes(int count)
{
    return sizeof(struct pool_bitmap) + ((count + 63u) / 64u) * sizeof(uint64_t);
}

static void bitmap_cache_add(struct swbwa_segment_tree *tree)
{
    unsigned char *cache = bitmap_caches[_MYID];
    size_t bytes = bitmap_bytes(tree->seg_size);
    if (!cache || !tree->bitmap_home || tree->bitmap != tree->bitmap_home ||
        bytes > POOL_BITMAP_CACHE_BYTES - bitmap_offsets[_MYID]) return;
    tree->bitmap = (struct pool_bitmap *)(cache + bitmap_offsets[_MYID]);
    memcpy(tree->bitmap, tree->bitmap_home, bytes);
    bitmap_offsets[_MYID] += bytes;
}

void swbwa_pool_bitmap_begin(unsigned mode)
{
    assert(!bitmap_caches[_MYID]);
    if (SWBWA_ALLOC_INITIALIZED && mode != bitmap_modes[_MYID])
        swbwa_cpe_fail(SWBWA_CPE_ERR_LDM_ALLOC_STATE, mode, bitmap_modes[_MYID], 104);
    bitmap_modes[_MYID] = mode;
    bitmap_offsets[_MYID] = 0;
    if (mode != 2) return;
    bitmap_caches[_MYID] = swbwa_ldm_alloc(POOL_BITMAP_CACHE_BYTES, 4);
    if (!SWBWA_ALLOC_INITIALIZED) return;
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; ++i)
        for (int j = 0; j < SWBWA_ALLOC_TREE_COUNTS[i]; ++j)
            bitmap_cache_add(SWBWA_ALLOC_TREES[i][j]);
}

void swbwa_pool_bitmap_end(void)
{
    unsigned char *cache = bitmap_caches[_MYID];
    if (!cache) return;
    if (SWBWA_ALLOC_INITIALIZED) {
        for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; ++i) {
            for (int j = 0; j < SWBWA_ALLOC_TREE_COUNTS[i]; ++j) {
                struct swbwa_segment_tree *t = SWBWA_ALLOC_TREES[i][j];
                if (t->bitmap == t->bitmap_home) continue;
                memcpy(t->bitmap_home, t->bitmap, bitmap_bytes(t->seg_size));
                t->bitmap = t->bitmap_home;
            }
        }
    }
    bitmap_caches[_MYID] = NULL;
    swbwa_ldm_release(cache, POOL_BITMAP_CACHE_BYTES);
}
#endif

#if SWBWA_LDM_UNIFIED
enum { POOL_CACHE_ENTRIES = 64 };
typedef struct {
    uintptr_t begin, end;
    unsigned short size_class, index;
} pool_cache_entry_t;
typedef struct {
    unsigned count;
    pool_cache_entry_t entry[POOL_CACHE_ENTRIES];
} pool_cache_t;
static pool_cache_t *pool_caches[SWBWA_CPE_COUNT];

static void pool_cache_add(struct swbwa_segment_tree *tree, int size_class, int index)
{
    pool_cache_t *cache = pool_caches[_MYID];
    unsigned i;
    if (!cache || cache->count == POOL_CACHE_ENTRIES) return;
    i = cache->count++;
    while (i && cache->entry[i - 1].begin > (uintptr_t)tree->start_address) {
        cache->entry[i] = cache->entry[i - 1];
        --i;
    }
    cache->entry[i].begin = (uintptr_t)tree->start_address;
    cache->entry[i].end = (uintptr_t)tree->end_address;
    cache->entry[i].size_class = size_class;
    cache->entry[i].index = index;
}

void swbwa_pool_cache_begin(void)
{
    pool_cache_t *cache;
    assert(!pool_caches[_MYID]);
    cache = swbwa_ldm_alloc(sizeof(*cache), 4);
    if (!cache) return;
    cache->count = 0;
    pool_caches[_MYID] = cache;
    if (!SWBWA_ALLOC_INITIALIZED) return;
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; ++i)
        for (int j = 0; j < SWBWA_ALLOC_TREE_COUNTS[i]; ++j)
            pool_cache_add(SWBWA_ALLOC_TREES[i][j], i, j);
}

void swbwa_pool_cache_end(void)
{
    pool_cache_t *cache = pool_caches[_MYID];
    pool_caches[_MYID] = NULL;
    if (cache) swbwa_ldm_release(cache, sizeof(*cache));
}
#endif

static int find_pool_tree(const void *ptr, int *size_class, int *tree_index)
{
    uintptr_t address = (uintptr_t)ptr;
#if SWBWA_LDM_UNIFIED
    /* Only descriptors are cached. Escaping SAM payload remains on the heap,
     * and a full/refused cache safely falls back to the original search. */
    const pool_cache_t *cache = pool_caches[_MYID];
    if (cache) {
        unsigned left = 0, right = cache->count;
        while (left < right) {
            unsigned mid = left + (right - left) / 2;
            if (cache->entry[mid].begin <= address) left = mid + 1;
            else right = mid;
        }
        if (left && address < cache->entry[left - 1].end) {
            *size_class = cache->entry[left - 1].size_class;
            *tree_index = cache->entry[left - 1].index;
            return 1;
        }
    }
#endif
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; ++i) {
        for (int j = 0; j < SWBWA_ALLOC_TREE_COUNTS[i]; ++j) {
            const struct swbwa_segment_tree *tree = SWBWA_ALLOC_TREES[i][j];
            if (address >= (uintptr_t)tree->start_address &&
                address < (uintptr_t)tree->end_address) {
                *size_class = i;
                *tree_index = j;
                return 1;
            }
        }
    }
    return 0;
}

/*
 * A single size class must not be able to swallow the whole pool. Tree sizes
 * double on every refill, so without a cap one class that briefly needs more
 * live blocks than its initial tree holds jumps straight to a multi-megabyte
 * allocation and starves everyone else.
 */
enum { SWBWA_ALLOC_MAX_TREE_BYTES = 1 << 20 };

long cpe_pool_high_water(void)
{
    return (long)pool_offsets[_MYID];
}

static long ldm_outstanding_bytes[SWBWA_CPE_COUNT];
static long ldm_peak_bytes[SWBWA_CPE_COUNT];
static long ldm_refusals[SWBWA_CPE_COUNT];

/*
 * The CPE stack shares LDM with this scratch, and overrunning it silently
 * overwrites saved return addresses instead of failing an allocation: the
 * core then returns to a garbage PC, which the hardware reports as an
 * unrelated integer/NPC overflow. ldm_malloc alone does not prevent that, so
 * cap the tracked scratch. This is a conservative budget, not a measurement
 * of total LDM usage: stack, static data and legacy raw allocations are extra.
 */

/*
 * Returns NULL when the request does not fit the budget or LDM refuses it;
 * every caller is expected to fall back to the heap. Occupancy varies with
 * the alignment work in flight, so a request that normally fits can
 * legitimately fail.
 */
void *swbwa_ldm_alloc(unsigned long bytes, int site)
{
    void *p;

#if SWBWA_CPE_LDM_MODE == SWBWA_CPE_LDM_OFF
    (void)bytes;
    (void)site;
    return NULL;
#elif !SWBWA_CPE_MANUAL_LDM
    /* Ablation: disable explicit scratch placement, not the new arena or
     * SIMD/algorithm changes. All old callers already have heap fallbacks. */
    if (site != SWBWA_LDM_AUTO_ARENA_SITE) {
#if SWBWA_LDM_UNIFIED
        return swbwa_ldm_scratch_try(bytes, (unsigned)site);
#else
        return NULL;
#endif
    }
#else
    (void)site;
#endif
    if (bytes > (unsigned long)SWBWA_LDM_SCRATCH_BUDGET_BYTES ||
        ldm_outstanding_bytes[_MYID] >
            SWBWA_LDM_SCRATCH_BUDGET_BYTES - (long)bytes) {
        ++ldm_refusals[_MYID];
        return NULL;
    }
    p = ldm_malloc(bytes);
    if (p == NULL) {
        ++ldm_refusals[_MYID];
        return NULL;
    }
    ldm_outstanding_bytes[_MYID] += (long)bytes;
    if (ldm_outstanding_bytes[_MYID] > ldm_peak_bytes[_MYID])
        ldm_peak_bytes[_MYID] = ldm_outstanding_bytes[_MYID];
    return p;
}

void swbwa_ldm_release(void *ptr, unsigned long bytes)
{
    if (ptr == NULL) return;
#if SWBWA_LDM_UNIFIED
    if (swbwa_ldm_scratch_release(ptr)) return;
#endif
    ldm_free(ptr, bytes);
    ldm_outstanding_bytes[_MYID] -= (long)bytes;
}

void swbwa_ldm_begin_batch(void)
{
    /* Preserve outstanding bytes so a leak is not hidden by the reset. */
    ldm_peak_bytes[_MYID] = ldm_outstanding_bytes[_MYID];
    ldm_refusals[_MYID] = 0;
}

long swbwa_ldm_outstanding(void)
{
    return ldm_outstanding_bytes[_MYID];
}

long swbwa_ldm_peak(void)
{
    return ldm_peak_bytes[_MYID];
}

long swbwa_ldm_refusals(void)
{
    return ldm_refusals[_MYID];
}

void set_big_buffer(char* buffer, long long t_size) {
    if (buffer == NULL || t_size <= 0) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_BAD_POOL, (long)t_size, 0, 0);
    }
    pool_starts[_MYID] = buffer + _MYID * t_size;
    pool_offsets[_MYID] = 0;
    pool_sizes[_MYID] = (size_t)t_size;
    allocator_initialized[_MYID] = 0;
    memset(tree_counts[_MYID], 0, sizeof(tree_counts[_MYID]));
    memset(next_tree_sizes[_MYID], 0, sizeof(next_tree_sizes[_MYID]));
}

void *l_calloc(size_t nmemb, size_t size) {
    size_t t_size;
    if (size != 0 && nmemb > SIZE_MAX / size) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_SIZE_OVERFLOW, (long)nmemb, (long)size, 0);
    }
    t_size = nmemb * size;
    if (t_size > pool_sizes[_MYID] - pool_offsets[_MYID]) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_EXHAUSTED, (long)t_size,
                       (long)(pool_sizes[_MYID] - pool_offsets[_MYID]),
                       (long)pool_offsets[_MYID]);
    }
    void *ptr = pool_starts[_MYID] + pool_offsets[_MYID];
    pool_offsets[_MYID] += t_size;
    memset(ptr, 0, t_size);
    return ptr;
}

static struct swbwa_segment_tree *build_segment_tree(int bind_length, int seg_size) {
    void *new_address = l_calloc(bind_length * seg_size, 1);
    struct swbwa_segment_tree *now_tree = (struct swbwa_segment_tree *) l_calloc(sizeof(struct swbwa_segment_tree), 1);
    now_tree->start_address = new_address;
    now_tree->end_address = now_tree->start_address + 1LL * bind_length * seg_size;
    now_tree->leaf_offset = seg_size - 1;
    now_tree->seg_size = seg_size;
    int tree_size = seg_size << 1;
    now_tree->tree = (int *) l_calloc(tree_size * sizeof(int), 1);
#if SWBWA_LDM_UNIFIED
    if (bitmap_modes[_MYID]) {
        /* Keep the payload and tree capacity unchanged. Only the free-index
         * representation changes; find_free still returns the lowest slot. */
        size_t bytes = bitmap_bytes(seg_size);
        now_tree->bitmap_home = bytes <= tree_size * sizeof(int)
            ? (struct pool_bitmap *)now_tree->tree : l_calloc(bytes, 1);
        now_tree->bitmap = now_tree->bitmap_home;
        if (seg_size % 64)
            now_tree->bitmap->words[seg_size / 64] =
                ~((UINT64_C(1) << (seg_size % 64)) - 1);
        bitmap_cache_add(now_tree);
        return now_tree;
    }
#endif
    for (int i = 1; i <= seg_size; ++i)
        now_tree->tree[i + now_tree->leaf_offset] = 0;
    for (int i = now_tree->leaf_offset; i; --i)
        now_tree->tree[i] = now_tree->tree[i << 1] + now_tree->tree[i << 1 | 1];
    return now_tree;
}

static void update_segment_tree(struct swbwa_segment_tree *now_tree, int pos, int value) {
#if SWBWA_LDM_UNIFIED
    if (now_tree->bitmap) {
        struct pool_bitmap *b = now_tree->bitmap;
        unsigned word = (unsigned)(pos - 1) / 64;
        uint64_t bit = UINT64_C(1) << ((pos - 1) % 64);
        assert(pos >= 1 && pos <= now_tree->seg_size);
        if (value == 1) {
            assert(!(b->words[word] & bit));
            b->words[word] |= bit;
            ++b->used;
        } else {
            assert(value == -1 && (b->words[word] & bit));
            b->words[word] &= ~bit;
            --b->used;
            if (word < b->first_word) b->first_word = word;
        }
        return;
    }
#endif
    for (int i = pos + now_tree->leaf_offset; i; i >>= 1)
        now_tree->tree[i] += value;
}

static int segment_allocation_state(struct swbwa_segment_tree *now_tree, int pos) {
#if SWBWA_LDM_UNIFIED
    if (now_tree->bitmap)
        return (now_tree->bitmap->words[(pos - 1) / 64] >> ((pos - 1) % 64)) & 1;
#endif
    return now_tree->tree[pos + now_tree->leaf_offset];
}

static int find_free_segment(struct swbwa_segment_tree *now_tree) {
#if SWBWA_LDM_UNIFIED
    if (now_tree->bitmap) {
        struct pool_bitmap *b = now_tree->bitmap;
        if (b->used == (unsigned)now_tree->seg_size) return -1;
        while (b->words[b->first_word] == UINT64_MAX) ++b->first_word;
        return (int)(b->first_word * 64 + __builtin_ctzll(~b->words[b->first_word]) + 1);
    }
#endif
    int l = 1, r = now_tree->seg_size;
    assert(now_tree->tree[1] <= now_tree->seg_size);
    if (now_tree->tree[1] == now_tree->seg_size) return -1;
    for (int i = 1; l != r;) {
        int mid = (l + r) / 2;
        if (now_tree->tree[i << 1] < (mid - l + 1)) {
            r = mid;
            i = i << 1;
        } else if (now_tree->tree[i << 1 | 1] < (r - mid)) {
            l = mid + 1;
            i = i << 1 | 1;
        } else {
            assert(0);
        }
    }
    return l;
}

static void allocator_init(void) {
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; i++) {
        SWBWA_ALLOC_TREE_COUNTS[i] = 0;
        struct swbwa_segment_tree *now_tree = build_segment_tree(block_sizes[i], initial_tree_sizes[i]);
        SWBWA_ALLOC_TREES[i][SWBWA_ALLOC_TREE_COUNTS[i]++] = now_tree;
#if SWBWA_LDM_UNIFIED
        pool_cache_add(now_tree, i, SWBWA_ALLOC_TREE_COUNTS[i] - 1);
#endif
        SWBWA_ALLOC_NEXT_TREE_SIZES[i] = initial_tree_sizes[i] << 1;
    }
    SWBWA_ALLOC_INITIALIZED = 1;
}

static void *segment_address(struct swbwa_segment_tree *now_tree, int bind_length, int seg_pos) {
    assert(seg_pos >= 1);
    return now_tree->start_address + bind_length * (seg_pos - 1);
}

static int segment_index(struct swbwa_segment_tree *now_tree, int bind_length, void *free_address) {
    return ((char *)free_address - now_tree->start_address) / bind_length + 1;
}

static void *pool_malloc(size_t size) {

    if (SWBWA_ALLOC_INITIALIZED == 0) allocator_init();

    /* size is clamped to [4, 2^18] by cpe_pool_malloc. */
    int length_type = 32 - __builtin_clz((unsigned int)(size - 1)) - 2;

    assert(length_type >= 0 && length_type < SWBWA_ALLOC_SIZE_CLASS_COUNT);

    int now_tree_num = SWBWA_ALLOC_TREE_COUNTS[length_type];
    int find_pos = -1;
    int first_zero_pos = -1;
    for (int i = 0; i < now_tree_num; i++) {
        struct swbwa_segment_tree *now_tree = SWBWA_ALLOC_TREES[length_type][i];
        first_zero_pos = find_free_segment(now_tree);
        if (first_zero_pos != -1) {
            find_pos = i;
            break;
        }
    }

    if (find_pos == -1) {
        if (SWBWA_ALLOC_TREE_COUNTS[length_type] >= SWBWA_ALLOC_MAX_TREES_PER_CLASS) {
            swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_TREE_LIMIT, length_type,
                           SWBWA_ALLOC_TREE_COUNTS[length_type], 0);
        }
        int next_size = SWBWA_ALLOC_NEXT_TREE_SIZES[length_type];

        /* Each segment also costs two ints of segment-tree index. */
        while (next_size > 1 &&
               (long long)next_size * (block_sizes[length_type] + 8) >
                   SWBWA_ALLOC_MAX_TREE_BYTES)
            next_size >>= 1;
        struct swbwa_segment_tree *now_tree = build_segment_tree(block_sizes[length_type], next_size);
        SWBWA_ALLOC_NEXT_TREE_SIZES[length_type] = next_size << 1;
        SWBWA_ALLOC_TREES[length_type][SWBWA_ALLOC_TREE_COUNTS[length_type]++] = now_tree;
#if SWBWA_LDM_UNIFIED
        pool_cache_add(now_tree, length_type, SWBWA_ALLOC_TREE_COUNTS[length_type] - 1);
#endif
        find_pos = SWBWA_ALLOC_TREE_COUNTS[length_type] - 1;
        first_zero_pos = 1;
    }

    assert(find_pos != -1 && first_zero_pos != -1);

    void *res_address = segment_address(SWBWA_ALLOC_TREES[length_type][find_pos], block_sizes[length_type], first_zero_pos);
    assert(segment_allocation_state(SWBWA_ALLOC_TREES[length_type][find_pos], first_zero_pos) == 0);
    update_segment_tree(SWBWA_ALLOC_TREES[length_type][find_pos], first_zero_pos, 1);

    return res_address;
}

void *cpe_pool_malloc(size_t size) {
    if(size < 4) size = 4;

    if (size > (1 << (SWBWA_ALLOC_SIZE_CLASS_COUNT + 1))) {
        return malloc(size);
    }
    return pool_malloc(size);
}

void cpe_pool_free(void *ptr) {
    if (ptr == NULL) return;
    if (SWBWA_ALLOC_INITIALIZED == 0) allocator_init();

    int pos1 = -1;
    int pos2 = -1;
    if (!find_pool_tree(ptr, &pos1, &pos2)) {
        free(ptr);
        return;
    }
    int free_seg_pos = segment_index(SWBWA_ALLOC_TREES[pos1][pos2], block_sizes[pos1], ptr);
    if (segment_allocation_state(SWBWA_ALLOC_TREES[pos1][pos2], free_seg_pos) != 1)
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_DOUBLE_FREE, pos1, pos2, free_seg_pos);
    update_segment_tree(SWBWA_ALLOC_TREES[pos1][pos2], free_seg_pos, -1);
}

void *cpe_pool_realloc(void *ptr, size_t size) {
    if (SWBWA_ALLOC_INITIALIZED == 0) allocator_init();
    if (ptr == NULL && size == 0) {
        return NULL;
    }

    if (ptr == NULL) {
        return cpe_pool_malloc(size);
    }

    if (size == 0) {
        cpe_pool_free(ptr);
        return NULL;
    }

    int pos1 = -1;
    int pos2 = -1;
    if (!find_pool_tree(ptr, &pos1, &pos2)) {
        return realloc(ptr, size);
    }

    size_t old_size = block_sizes[pos1];
    if (size > old_size) {
        void *new_ptr = cpe_pool_malloc(size);
        if (new_ptr == NULL) return NULL;
        memcpy(new_ptr, ptr, old_size);
        int free_seg_pos = segment_index(SWBWA_ALLOC_TREES[pos1][pos2], block_sizes[pos1], ptr);
        assert(segment_allocation_state(SWBWA_ALLOC_TREES[pos1][pos2], free_seg_pos) == 1);
        update_segment_tree(SWBWA_ALLOC_TREES[pos1][pos2], free_seg_pos, -1);
        return new_ptr;
    }
    return ptr;
}

char *cpe_pool_strdup(const char *s) {
    char *copy = cpe_pool_malloc(strlen(s) + 1);
    if (copy) {
        strcpy(copy, s);
    }
    return copy;
}

void *wrap_calloc(size_t nmemb, size_t size,
                  const char *file, unsigned int line, const char *func)
{
    size_t bytes;
    void *p;

    if (size != 0 && nmemb > SIZE_MAX / size) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_POOL_SIZE_OVERFLOW, (long)nmemb, (long)size, line);
    }
    bytes = nmemb * size;
    p = cpe_pool_malloc(bytes);
    if (bytes > 0 && p == NULL) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)bytes, line, 0);
    }
    memset(p, 0, bytes);
    return p;
}

void *wrap_malloc(size_t size,
                  const char *file, unsigned int line, const char *func)
{
    void *p = cpe_pool_malloc(size);
    if (size > 0 && p == NULL) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)size, line, 0);
    }
    return p;
}

void *wrap_realloc(void *ptr, size_t size,
                   const char *file, unsigned int line, const char *func)
{
    void *p = cpe_pool_realloc(ptr, size);
    if (size > 0 && p == NULL) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)size, line, 0);
    }
    return p;
}

char *wrap_strdup(const char *s,
                  const char *file, unsigned int line, const char *func)
{
    char *p = cpe_pool_strdup(s);
    if (p == NULL) {
        swbwa_cpe_fail(SWBWA_CPE_ERR_ALLOC_FAILED, (long)strlen(s), line, 0);
    }
    return p;
}

void wrap_free(void *ptr, const char *file, unsigned int line, const char *func)
{
    (void)file;
    (void)line;
    (void)func;
    if (ptr != NULL) cpe_pool_free(ptr);
}
