/* Host-only tests of the actual allocator, with a stubbed SDK LDM allocator. */
#include <assert.h>
#include <limits.h>
#include <setjmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swbwa_config.h"
#include "malloc_wrap.h"

int swbwa_test_cpe_id;
static int refuse_ldm, expected_error;
static jmp_buf failure;

void *ldm_malloc(unsigned long bytes)
{
    return refuse_ldm ? NULL : malloc(bytes);
}

void ldm_free(void *ptr, unsigned long bytes)
{
    (void)bytes;
    free(ptr);
}

void swbwa_cpe_fail(int code, long a, long b, long c)
{
    (void)a; (void)b; (void)c;
    assert(code == expected_error);
    longjmp(failure, 1);
}

void *l_calloc(size_t nmemb, size_t size);
void *cpe_pool_malloc(size_t size);
void cpe_pool_free(void *ptr);
void *cpe_pool_realloc(void *ptr, size_t size);

int main(void)
{
    void *p, *q;
    char *pool = malloc(16 << 20);
    void *blocks[32];
    int i;

    assert(pool != NULL);
    swbwa_ldm_begin_batch();
    p = swbwa_ldm_alloc(40 << 10, 0);
    assert(p != NULL && swbwa_ldm_outstanding() == 40 << 10);
    assert(swbwa_ldm_alloc(1, 0) == NULL);
    assert(swbwa_ldm_alloc(ULONG_MAX, 0) == NULL);
    assert(swbwa_ldm_refusals() == 2);
    swbwa_ldm_begin_batch();
    assert(swbwa_ldm_peak() == 40 << 10 && swbwa_ldm_refusals() == 0);
    swbwa_ldm_release(p, 40 << 10);
    swbwa_ldm_begin_batch();
    assert(swbwa_ldm_outstanding() == 0 && swbwa_ldm_peak() == 0);
    refuse_ldm = 1;
    assert(swbwa_ldm_alloc(8, 0) == NULL && swbwa_ldm_outstanding() == 0);
    assert(swbwa_ldm_refusals() == 1);
    refuse_ldm = 0;
    swbwa_test_cpe_id = 1;
    assert(swbwa_ldm_refusals() == 0 && swbwa_ldm_outstanding() == 0);
    swbwa_test_cpe_id = 0;

    set_big_buffer(pool, 16 << 20);
    p = cpe_pool_malloc(8);
    memset(p, 0x5a, 8);
    q = cpe_pool_realloc(p, 512);
    assert(q != NULL && memcmp(q, "ZZZZZZZZ", 8) == 0);
    cpe_pool_free(q);
    /* Exercise the bounded refills of the largest size class. */
    for (i = 0; i < 32; ++i) {
        long before = cpe_pool_high_water();
        blocks[i] = cpe_pool_malloc(1 << 18);
        assert(cpe_pool_high_water() - before < (1 << 20) + 128);
        memset(blocks[i], i, 1 << 18);
    }
    for (i = 0; i < 32; ++i) {
        assert(((unsigned char *)blocks[i])[0] == i);
        cpe_pool_free(blocks[i]);
    }
    expected_error = SWBWA_CPE_ERR_POOL_DOUBLE_FREE;
    if (setjmp(failure) == 0) {
        cpe_pool_free(blocks[0]);
        assert(0);
    }
    expected_error = SWBWA_CPE_ERR_POOL_SIZE_OVERFLOW;
    if (setjmp(failure) == 0) {
        l_calloc(SIZE_MAX, 2);
        assert(0);
    }
    expected_error = SWBWA_CPE_ERR_POOL_EXHAUSTED;
    if (setjmp(failure) == 0) {
        l_calloc(SIZE_MAX, 1);
        assert(0);
    }
    free(pool);
    puts("CPE allocator host tests passed");
    return 0;
}
