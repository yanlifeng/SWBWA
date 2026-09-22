#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <sys/stat.h>
#include <unistd.h>

#include <slave.h>

#if SWBWA_ENABLE_CPE_MALLOC_WRAPPER
/* Don't wrap ourselves */
#  undef SWBWA_ENABLE_CPE_MALLOC_WRAPPER
#endif

#include "malloc_wrap.h"

#define SWBWA_ALLOC_GLOBAL

struct swbwa_segment_tree {
    int *tree;
    int seg_size;
    int leaf_offset;
    char *start_address;
    char *end_address;
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

SWBWA_ALLOC_GLOBAL int allocator_initialized[SWBWA_CPE_COUNT] = {0};
SWBWA_ALLOC_GLOBAL struct swbwa_segment_tree
    *allocator_trees[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT]
                    [SWBWA_ALLOC_MAX_TREES_PER_CLASS];
SWBWA_ALLOC_GLOBAL int tree_counts[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT];
SWBWA_ALLOC_GLOBAL int next_tree_sizes[SWBWA_CPE_COUNT][SWBWA_ALLOC_SIZE_CLASS_COUNT];

#define SWBWA_ALLOC_INITIALIZED allocator_initialized[_MYID]
#define SWBWA_ALLOC_TREES allocator_trees[_MYID]
#define SWBWA_ALLOC_TREE_COUNTS tree_counts[_MYID]
#define SWBWA_ALLOC_NEXT_TREE_SIZES next_tree_sizes[_MYID]

static char *pool_starts[SWBWA_CPE_COUNT];
static size_t pool_offsets[SWBWA_CPE_COUNT];
static size_t pool_sizes[SWBWA_CPE_COUNT];

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
    if (site != SWBWA_LDM_AUTO_ARENA_SITE) return NULL;
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
    for (int i = 1; i <= seg_size; ++i)
        now_tree->tree[i + now_tree->leaf_offset] = 0;
    for (int i = now_tree->leaf_offset; i; --i)
        now_tree->tree[i] = now_tree->tree[i << 1] + now_tree->tree[i << 1 | 1];
    return now_tree;
}

static void update_segment_tree(struct swbwa_segment_tree *now_tree, int pos, int value) {
    for (int i = pos + now_tree->leaf_offset; i; i >>= 1)
        now_tree->tree[i] += value;
}

static int segment_allocation_state(struct swbwa_segment_tree *now_tree, int pos) {
    return now_tree->tree[pos + now_tree->leaf_offset];
}

static int find_free_segment(struct swbwa_segment_tree *now_tree) {
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

    int length_type = (int) ceil(log2(size)) - 2;

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
    if (SWBWA_ALLOC_INITIALIZED == 0) allocator_init();

    int pos1 = -1;
    int pos2 = -1;
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; i++) {
        int now_tree_num = SWBWA_ALLOC_TREE_COUNTS[i];
        for (int j = 0; j < now_tree_num; j++) {
            uintptr_t p1 = (uintptr_t)SWBWA_ALLOC_TREES[i][j]->start_address;
            uintptr_t p2 = (uintptr_t)SWBWA_ALLOC_TREES[i][j]->end_address;
            uintptr_t p3 = (uintptr_t)ptr;
            if (p3 >= p1 && p3 < p2) {
                pos1 = i;
                pos2 = j;
                break;
            }
        }
        if (pos1 != -1) break;
    }
    if (pos1 == -1 && pos2 == -1) {
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
    for (int i = 0; i < SWBWA_ALLOC_SIZE_CLASS_COUNT; i++) {
        int now_tree_num = SWBWA_ALLOC_TREE_COUNTS[i];
        for (int j = 0; j < now_tree_num; j++) {
            uintptr_t p1 = (uintptr_t)SWBWA_ALLOC_TREES[i][j]->start_address;
            uintptr_t p2 = (uintptr_t)SWBWA_ALLOC_TREES[i][j]->end_address;
            uintptr_t p3 = (uintptr_t)ptr;
            if (p3 >= p1 && p3 < p2) {
                pos1 = i;
                pos2 = j;
                break;
            }
        }
        if (pos1 != -1) break;
    }

    if (pos1 == -1 && pos2 == -1) {
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

void swbwa_cpe_malloc_stats_init(void) {
    char folder[256] = "cpe_malloc_stats";
    if (access(folder, F_OK) == -1) {
        mkdir(folder, 0777);
    }
    //timer_flag = 1;
    if (SWBWA_ALLOC_INITIALIZED == 0) allocator_init();
}

long cal(int id, int n) {
    long res = 0;
    long now = initial_tree_sizes[id];
    for (int i = 0; i < n; i++) {
        res += now;
        now = now << 1;
    }
    return res;
}

void swbwa_cpe_malloc_stats_print(void) {

    int my_cpe_id = _MYID;
    char filename[256];
    FILE *fp;

    sprintf(filename, "cpe_malloc_stats/cpe_malloc_rank%06d.dat", my_cpe_id);
    fp = fopen(filename, "w");

    fprintf(fp, "list size info : \n");
    fprintf(fp, "[0B ,4B]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[0], cal(0, SWBWA_ALLOC_TREE_COUNTS[0]));
    fprintf(fp, "(4B ,8B]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[1], cal(1, SWBWA_ALLOC_TREE_COUNTS[1]));
    fprintf(fp, "(8B ,16B]      == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[2], cal(2, SWBWA_ALLOC_TREE_COUNTS[2]));
    fprintf(fp, "(16B ,32B]     == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[3], cal(3, SWBWA_ALLOC_TREE_COUNTS[3]));
    fprintf(fp, "(32B ,64B]     == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[4], cal(4, SWBWA_ALLOC_TREE_COUNTS[4]));
    fprintf(fp, "(64B ,128B]    == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[5], cal(5, SWBWA_ALLOC_TREE_COUNTS[5]));
    fprintf(fp, "(128B ,256B]   == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[6], cal(6, SWBWA_ALLOC_TREE_COUNTS[6]));
    fprintf(fp, "(256B ,512B]   == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[7], cal(7, SWBWA_ALLOC_TREE_COUNTS[7]));
    fprintf(fp, "(512B ,1K]     == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[8], cal(8, SWBWA_ALLOC_TREE_COUNTS[8]));
    fprintf(fp, "(1K ,2K]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[9], cal(9, SWBWA_ALLOC_TREE_COUNTS[9]));
    fprintf(fp, "(2K ,4K]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[10], cal(10, SWBWA_ALLOC_TREE_COUNTS[10]));
    fprintf(fp, "(4K ,8K]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[11], cal(11, SWBWA_ALLOC_TREE_COUNTS[11]));
    fprintf(fp, "(8K ,16K]      == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[12], cal(12, SWBWA_ALLOC_TREE_COUNTS[12]));
    fprintf(fp, "(16K ,32K]     == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[13], cal(13, SWBWA_ALLOC_TREE_COUNTS[13]));
    fprintf(fp, "(32K ,64K]     == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[14], cal(14, SWBWA_ALLOC_TREE_COUNTS[14]));
    fprintf(fp, "(64K ,128K]    == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[15], cal(15, SWBWA_ALLOC_TREE_COUNTS[15]));
    //fprintf(fp, "(64K ,-]       == list size %8d  tot %12ld ele\n", SWBWA_ALLOC_TREE_COUNTS[13], cal(SWBWA_ALLOC_TREE_COUNTS[13]));

    //fprintf(fp, "align address infp : \n");

    //for(int i = 0; i < test_address_size; i++) {
    //    fprintf(fp, "%p   ",valloc_addss[i]);
    //}
    //fprintf(fp, "\n\n");

    //fprintf(fp, "align address infp : \n");
    //for(int i = 0; i < test_address_size; i++) {
    //    fprintf(fp, "%p   ", memalign_addss[i]);
    //}
    //fprintf(fp, "\n\n");



    fclose(fp);
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
