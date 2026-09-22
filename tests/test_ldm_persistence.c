/* Build the same source for MPE and CPE. No alignment or input data required. */
#include <stdint.h>
#include <stddef.h>

#ifndef TEST_LDM_BYTES
#define TEST_LDM_BYTES 32768
#endif
#ifndef TEST_LDM_STACK_BYTES
#define TEST_LDM_STACK_BYTES 8192
#endif
#ifndef TEST_LDM_MIN_STACK_GAP
#define TEST_LDM_MIN_STACK_GAP 0
#endif
enum { CORES = 384, BYTES = TEST_LDM_BYTES, PHASES = 32 };
struct record {
    uintptr_t address, initial_address, stack_address;
    unsigned checks, failures, overlaps, tls_changes;
};
struct task { unsigned phase; struct record record[CORES]; };

#ifdef TEST_LDM_CPE
#include <slave.h>
#include <string.h>
static __thread unsigned saved_phase;

static unsigned char pattern(unsigned core, unsigned phase, unsigned i)
{
    return (unsigned char)(core * 17u + phase * 31u + i * 7u + (i >> 8));
}

static __attribute__((noinline)) void use_stack(unsigned value)
{
    volatile unsigned char scratch[TEST_LDM_STACK_BYTES];
    for (unsigned i = 0; i < sizeof(scratch); ++i)
        scratch[i] = (unsigned char)(value + i);
    asm volatile("" : : "r"(&scratch[0]) : "memory");
}

void test_ldm_persistence(struct task *task)
{
    unsigned id = (unsigned)_MYID, phase = task->phase;
    struct record *r = &task->record[id];
    unsigned char *p = (unsigned char *)r->address;
    if (phase == 0) {
        p = ldm_malloc(BYTES);
        r->address = (uintptr_t)p;
        r->initial_address = (uintptr_t)p;
        r->stack_address = (uintptr_t)&p;
        if (!p) { ++r->failures; goto done; }
#if TEST_LDM_MIN_STACK_GAP > 0
        if ((uintptr_t)p >= r->stack_address ||
            BYTES > r->stack_address - (uintptr_t)p ||
            r->stack_address - (uintptr_t)p - BYTES < TEST_LDM_MIN_STACK_GAP) {
            ldm_free(p, BYTES);
            r->address = 0;
            ++r->failures;
            goto done;
        }
#endif
    } else {
        if (!p || r->overlaps) goto done;
        if (saved_phase != phase - 1) ++r->tls_changes;
        for (unsigned i = 0; i < BYTES; ++i)
            if (p[i] != pattern(id, phase - 1, i)) { ++r->failures; break; }
        ++r->checks;
    }
    use_stack(phase);
    {
        unsigned char *extra = ldm_malloc(4096);
        if (!extra) { ++r->failures; goto done; }
        if ((uintptr_t)extra < (uintptr_t)p + BYTES &&
            (uintptr_t)p < (uintptr_t)extra + 4096) {
            ++r->overlaps;
            goto done;
        }
        memset(extra, (int)phase, 4096);
        ldm_free(extra, 4096);
    }
    if (phase != 0)
        for (unsigned i = 0; i < BYTES; ++i)
            if (p[i] != pattern(id, phase - 1, i)) { ++r->failures; break; }
    if (phase == PHASES) {
        ldm_free(p, BYTES);
        r->address = 0;
    } else {
        for (unsigned i = 0; i < BYTES; ++i) p[i] = pattern(id, phase, i);
    }
    saved_phase = phase;
done:
    asm volatile("memb\n\t" ::: "memory");
    flush_slave_cache();
}
#else
#include <athread.h>
#include <stdio.h>
#include <stdlib.h>

extern void SLAVE_FUN(test_ldm_persistence)();
static __uncached struct task task;
int main(void)
{
    unsigned failures = 0, overlaps = 0, checks = 0, tls_changes = 0;
    uintptr_t min_gap = UINTPTR_MAX;
    athread_init_cgs();
    for (unsigned phase = 0; phase <= PHASES; ++phase) {
        task.phase = phase;
        __real_athread_spawn_cgs((void *)slave_test_ldm_persistence, &task, 1);
        athread_join_cgs();
    }
    for (unsigned i = 0; i < CORES; ++i) {
        struct record *r = &task.record[i];
        failures += r->failures; overlaps += r->overlaps;
        checks += r->checks; tls_changes += r->tls_changes;
        if (r->address) ++failures;
        if (r->initial_address && r->stack_address >= r->initial_address + BYTES) {
            uintptr_t gap = r->stack_address - r->initial_address - BYTES;
            if (gap < min_gap) min_gap = gap;
        }
        if (i == 0)
            printf("LDM_LAYOUT initial=%#lx stack=%#lx request=%u stack_test=%u\n",
                   (unsigned long)r->initial_address, (unsigned long)r->stack_address,
                   BYTES, TEST_LDM_STACK_BYTES);
        if (r->failures || r->overlaps || r->checks != PHASES || r->address)
            printf("core=%u checks=%u failures=%u overlaps=%u tls_changes=%u address=%#lx\n",
                   i, r->checks, r->failures, r->overlaps, r->tls_changes,
                   (unsigned long)r->address);
    }
    printf("LDM_PERSIST cores=%d bytes=%d spawns=%d checks=%u failures=%u overlaps=%u tls_changes=%u\n",
           CORES, BYTES, PHASES + 1, checks, failures, overlaps, tls_changes);
    printf("LDM_STACK min_initial_gap=%lu required=%u\n",
           (unsigned long)min_gap, TEST_LDM_MIN_STACK_GAP);
    return failures || overlaps || checks != CORES * PHASES ? EXIT_FAILURE : EXIT_SUCCESS;
}
#endif
