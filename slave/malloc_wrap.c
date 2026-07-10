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

void set_big_buffer(char* buffer, long long t_size) {
    if (buffer == NULL || t_size <= 0) {
        printf("invalid CPE allocator pool: buffer=%p size=%lld\n", buffer, t_size);
        exit(EXIT_FAILURE);
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
        printf("CPE allocator size overflow: %zu * %zu\n", nmemb, size);
        exit(EXIT_FAILURE);
    }
    t_size = nmemb * size;
    if (pool_offsets[_MYID] + t_size > pool_sizes[_MYID]) {
        printf("CPE allocator pool exhausted: requested=%zu available=%zu\n",
               t_size, pool_sizes[_MYID] - pool_offsets[_MYID]);
        exit(EXIT_FAILURE);
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
            printf("CPE allocator tree limit reached for size class %d\n", length_type);
            exit(EXIT_FAILURE);
        }
        struct swbwa_segment_tree *now_tree = build_segment_tree(block_sizes[length_type], SWBWA_ALLOC_NEXT_TREE_SIZES[length_type]);
        SWBWA_ALLOC_NEXT_TREE_SIZES[length_type] = SWBWA_ALLOC_NEXT_TREE_SIZES[length_type] << 1;
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
    assert(segment_allocation_state(SWBWA_ALLOC_TREES[pos1][pos2], free_seg_pos) == 1);
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
        printf(
                "[%s] calloc size overflow at %s line %u\n",
                func, file, line);
        exit(EXIT_FAILURE);
    }
    bytes = nmemb * size;
    p = cpe_pool_malloc(bytes);
    if (bytes > 0 && p == NULL) {
        printf("[%s] failed to calloc %zu bytes at %s line %u: %s\n",
               func, bytes, file, line, strerror(errno));
        exit(EXIT_FAILURE);
    }
    memset(p, 0, bytes);
    return p;
}

void *wrap_malloc(size_t size,
                  const char *file, unsigned int line, const char *func)
{
    void *p = cpe_pool_malloc(size);
    if (size > 0 && p == NULL) {
        printf("[%s] failed to malloc %zu bytes at %s line %u: %s\n",
               func, size, file, line, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return p;
}

void *wrap_realloc(void *ptr, size_t size,
                   const char *file, unsigned int line, const char *func)
{
    void *p = cpe_pool_realloc(ptr, size);
    if (size > 0 && p == NULL) {
        printf("[%s] failed to realloc %zu bytes at %s line %u: %s\n",
               func, size, file, line, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return p;
}

char *wrap_strdup(const char *s,
                  const char *file, unsigned int line, const char *func)
{
    char *p = cpe_pool_strdup(s);
    if (p == NULL) {
        printf("[%s] failed to strdup %zu bytes at %s line %u: %s\n",
               func, strlen(s), file, line, strerror(errno));
        exit(EXIT_FAILURE);
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
