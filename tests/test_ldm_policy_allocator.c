#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SWBWA_ALLOC_IMPLEMENTATION
#include "malloc_wrap.h"

int swbwa_test_cpe_id;
static void *backing;
static size_t backing_bytes;
static int refuse;
static uint32_t state = 371829;

void *ldm_malloc(unsigned long bytes)
{
    assert(!backing);
    if (refuse) return NULL;
    backing = malloc(bytes);
    backing_bytes = bytes;
    return backing;
}

void ldm_free(void *ptr, unsigned long bytes)
{
    assert(ptr == backing && bytes == backing_bytes);
    free(ptr);
    backing = NULL;
}

void swbwa_cpe_fail(int code, long a, long b, long c)
{
    fprintf(stderr, "CPE error %d %ld %ld %ld\n", code, a, b, c);
    abort();
}

static uint32_t next_random(void)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static int local(const void *ptr)
{
    uintptr_t p = (uintptr_t)ptr, base = (uintptr_t)backing;
    return backing && p >= base && p - base < backing_bytes;
}

static void parse(const char *text, unsigned *values)
{
    for (int i = 0; i < SWBWA_LDM_SITE_COUNT; ++i) {
        char *end;
        unsigned long value = strtoul(text, &end, 10);
        assert(end != text && value <= SWBWA_CPE_LDM_BYTES);
        values[i] = (unsigned)value;
        assert(*end == (i + 1 == SWBWA_LDM_SITE_COUNT ? '\0' : ','));
        text = end + (i + 1 < SWBWA_LDM_SITE_COUNT);
    }
}

static void check(const unsigned char *p, size_t size, unsigned char value)
{
    for (size_t i = 0; i < size; ++i) assert(p[i] == value);
}

int main(int argc, char **argv)
{
    enum { OBJECTS = 96, UPDATES = 20000 };
    swbwa_ldm_policy_t policy = {0};
    unsigned char *objects[OBJECTS] = {0};
    size_t sizes[OBJECTS] = {0};
    unsigned sites[OBJECTS] = {0};
    char *heap = malloc(128u << 20);
    assert((argc == 1 || argc == 3) && heap);
    set_big_buffer(heap, 128u << 20);
    swbwa_ldm_policy_default(&policy);
    if (argc == 3) {
        parse(argv[1], policy.cap);
        parse(argv[2], policy.reserve);
    }
    for (refuse = 0; refuse < 2; ++refuse) {
        for (unsigned profile = 0; profile < 2; ++profile) {
            policy.profile = profile;
            swbwa_ldm_set_policy(&policy);
            swbwa_ldm_allocator_begin();
            void *smem = swbwa_ldm_alloc(16384, 1);
            if (smem) memset(smem, 0x73, 16384);
            for (int k = 0; k < UPDATES; ++k) {
                unsigned index = next_random() % OBJECTS;
                unsigned char byte = (unsigned char)(index + 1);
                if (objects[index]) {
                    check(objects[index], sizes[index], byte);
                    if (next_random() % 3 == 0) {
                        size_t bytes = next_random() % 24000 + 1;
                        objects[index] = swbwa_auto_realloc(objects[index], bytes, sites[index]);
                        check(objects[index], sizes[index] < bytes ? sizes[index] : bytes, byte);
                        memset(objects[index], byte, bytes);
                        sizes[index] = bytes;
                    } else {
                        swbwa_auto_free(objects[index]);
                        objects[index] = NULL;
                    }
                } else {
                    unsigned site = next_random() % 5;
                    size_t bytes = next_random() % 20000 + 1;
                    objects[index] = swbwa_auto_malloc(bytes, site);
                    assert(objects[index]);
                    if (site == 0 || bytes > policy.cap[site] || refuse)
                        assert(!local(objects[index]));
                    memset(objects[index], byte, bytes);
                    sizes[index] = bytes;
                    sites[index] = site;
                }
            }
            for (unsigned i = 0; i < OBJECTS; ++i) {
                if (objects[i]) {
                    check(objects[i], sizes[i], (unsigned char)(i + 1));
                    swbwa_auto_free(objects[i]);
                    objects[i] = NULL;
                }
            }
            if (smem) {
                check(smem, 16384, 0x73);
                swbwa_ldm_release(smem, 16384);
            }
            unsigned char *escaped = swbwa_auto_malloc(65536, SWBWA_LDM_SITE_OTHER);
            assert(!local(escaped));
            memset(escaped, 0x42, 65536);
            swbwa_ldm_allocator_suspend();
            swbwa_ldm_allocator_resume();
            check(escaped, 65536, 0x42);
            swbwa_auto_free(escaped);
            swbwa_ldm_allocator_end();
            assert(!backing && swbwa_ldm_outstanding() == 0);
        }
    }
    free(heap);
    puts("PASS candidate: contents/realloc/spill/fragmentation/hints/refusal/lifetime");
    return 0;
}
