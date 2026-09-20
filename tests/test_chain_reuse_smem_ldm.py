"""Exercise SMEM/chain overlay ownership using actual collection and BWT code."""

import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HARNESS = r'''
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "malloc_wrap.h"

static size_t allocations, releases;
#if SWBWA_TEST_POOL_ALLOC
static void *count_malloc(size_t n)
{ ++allocations; return wrap_malloc(n, __FILE__, __LINE__, __func__); }
static void *count_calloc(size_t n, size_t s)
{ ++allocations; return wrap_calloc(n, s, __FILE__, __LINE__, __func__); }
static void *count_realloc(void *p, size_t n)
{ ++allocations; return wrap_realloc(p, n, __FILE__, __LINE__, __func__); }
static void count_free(void *p)
{ if (p) ++releases; wrap_free(p, __FILE__, __LINE__, __func__); }
#else
static void *count_malloc(size_t n) { ++allocations; return malloc(n); }
static void *count_calloc(size_t n, size_t s)
{ ++allocations; return calloc(n, s); }
static void *count_realloc(void *p, size_t n)
{ ++allocations; return realloc(p, n); }
static void count_free(void *p) { if (p) ++releases; free(p); }
#endif

#define malloc count_malloc
#define calloc count_calloc
#define realloc count_realloc
#define free count_free
#include "src/slave/bwamem.c"
#include "src/slave/bwt.c"
#undef malloc
#undef calloc
#undef realloc
#undef free

int swbwa_test_cpe_id;
static int refuse_ldm;
void *ldm_malloc(unsigned long n) { return refuse_ldm ? NULL : malloc(n); }
void ldm_free(void *p, unsigned long n) { (void)n; free(p); }
void swbwa_cpe_fail(int code, long a, long b, long c)
{
    fprintf(stderr, "unexpected CPE error: %d %ld %ld %ld\n", code, a, b, c);
    abort();
}

static void test_growth(smem_aux_t *aux)
{
    mem_opt_t opt = {0};
    mem_chain_t chain = {0};
    mem_seed_t seed = {0};
    int i;
    opt.w = 100;
    opt.max_chain_gap = 10000;
    aux->chain_arena_used = 0;
    chain.n = 1;
    chain.m = 4;
    chain.rid = 7;
    chain.seeds = smem_chain_alloc(aux, chain.m * sizeof(mem_seed_t), 1);
    seed.rbeg = 1000; seed.qbeg = 0; seed.len = seed.score = 19;
    chain.seeds[0] = seed;
    for (i = 1; i < 1024; ++i) {
        seed.rbeg = 1000 + i * 10; seed.qbeg = i * 10;
        assert(test_and_merge(&opt, 1000000, aux, &chain, &seed, 7));
        assert(chain.n == i + 1);
    }
    for (i = 0; i < chain.n; ++i) {
        assert(chain.seeds[i].rbeg == 1000 + i * 10);
        assert(chain.seeds[i].qbeg == i * 10);
        assert(chain.seeds[i].len == 19 && chain.seeds[i].score == 19);
    }
    /* Containment and reference rejection must not append a seed. */
    assert(test_and_merge(&opt, 1000000, aux, &chain, &seed, 7));
    assert(!test_and_merge(&opt, 1000000, aux, &chain, &seed, 8));
    assert(chain.n == 1024);
    smem_chain_free(aux, chain.seeds);
}

static void test_boundaries(smem_aux_t *aux)
{
    unsigned char *p, *q;
    size_t before;
    size_t capacity = aux->chain_arena ? smem_chain_arena_capacity(aux) :
                                       SWBWA_CHAIN_LDM_ARENA_BYTES;
    aux->chain_arena_used = 0;
    p = smem_chain_alloc(aux, capacity, 1);
    for (size_t i = 0; i < capacity; ++i) assert(p[i] == 0);
    if (aux->chain_arena) {
        assert(smem_chain_arena_owns(aux, p));
        assert(!smem_chain_arena_owns(aux, p + capacity));
    }
    q = smem_chain_alloc(aux, 8, 0);
    assert(!smem_chain_arena_owns(aux, q));
    memset(q, 0xa5, 8);
    before = releases;
    smem_chain_free(aux, q);
    assert(releases == before + 1);
    smem_chain_free(aux, p);
    aux->chain_arena_used = 0;
    p = smem_chain_alloc(aux, 1, 0);
    q = smem_chain_alloc(aux, 1, 0);
    if (aux->chain_arena) {
        assert(((uintptr_t)p & 7) == 0 && ((uintptr_t)q & 7) == 0);
        assert(q == p + 8);
    }
    smem_chain_free(aux, p);
    smem_chain_free(aux, q);
}

static size_t smem_bytes(void)
{
    return 2 * (size_t)SWBWA_SMEM_LDM_TMP_CAPACITY * sizeof(bwtintv_t);
}

static int is_overlay(const smem_aux_t *aux)
{
    return aux->chain_arena &&
           aux->chain_arena == (unsigned char *)aux->ldm_tmp_storage;
}

static void check_ownership(smem_aux_t *aux, long context_bytes)
{
    int expected_overlay = SWBWA_CHAIN_REUSE_SMEM_LDM && aux->ldm_tmp_storage;
    long expected = context_bytes;
    assert(is_overlay(aux) == expected_overlay);
    if (is_overlay(aux)) {
        size_t expected_capacity = smem_bytes() < SWBWA_CHAIN_LDM_ARENA_BYTES ?
                                  smem_bytes() : SWBWA_CHAIN_LDM_ARENA_BYTES;
        assert(smem_chain_arena_capacity(aux) == expected_capacity);
    } else if (aux->chain_arena) {
        assert(smem_chain_arena_capacity(aux) == SWBWA_CHAIN_LDM_ARENA_BYTES);
    }
    if (aux->ldm_tmp_storage) expected += smem_bytes();
    if (aux->chain_arena && !is_overlay(aux))
        expected += SWBWA_CHAIN_LDM_ARENA_BYTES;
    assert(swbwa_ldm_outstanding() == expected);
    assert(swbwa_ldm_peak() <= (40 << 10));
}

/* Exact packed BWT of A^512 T^512. Its reverse complement is itself. */
static void init_test_bwt(bwt_t *bwt)
{
    bwtint_t counts[4] = {0};
    memset(bwt, 0, sizeof(*bwt));
    bwt->primary = 1;
    bwt->seq_len = 1024;
    bwt->L2[1] = bwt->L2[2] = bwt->L2[3] = 512;
    bwt->L2[4] = 1024;
    bwt->bwt_size = 8 * 16;
    bwt->bwt = calloc(bwt->bwt_size, sizeof(uint32_t));
    assert(bwt->bwt);
    for (int i = 0; i < 1024; ++i) {
        uint32_t *block = bwt->bwt + (i / 128) * 16;
        int base = (i == 0 || (i >= 512 && i < 1023)) ? 3 : 0;
        if (i % 128 == 0) memcpy(block, counts, sizeof(counts));
        block[8 + (i % 128) / 16] |= (uint32_t)base << ((15 - i % 16) * 2);
        ++counts[base];
    }
    bwt_gen_cnt_table(bwt);
}

static void equal_intervals(const bwtintv_v *a, const bwtintv_v *b)
{
    assert(a->n == b->n);
    for (size_t i = 0; i < a->n; ++i) {
        assert(a->a[i].info == b->a[i].info);
        /* Backward extension leaves x[1] unspecified. Chaining consumes
         * the forward interval x[0], its size x[2], and query coordinates. */
        assert(a->a[i].x[0] == b->a[i].x[0]);
        assert(a->a[i].x[2] == b->a[i].x[2]);
    }
}

static void test_collection_phases(smem_aux_t *aux)
{
    const int lengths[] = {19, 63, 64, 65, 75, 150, 255, 256, 383, 384, 511, 512};
    bwt_t bwt;
    mem_opt_t opt = {0};
    smem_aux_t *reference = smem_aux_init(0);
    uint8_t query[512];
    int nonempty = 0;
    opt.min_seed_len = 19;
    opt.split_factor = 1.5f;
    opt.split_width = 10;
    opt.max_mem_intv = 20;
    init_test_bwt(&bwt);
    for (int round = 0; round < 3; ++round) {
        for (size_t j = 0; j < sizeof(lengths) / sizeof(lengths[0]); ++j) {
            int len = lengths[j];
            void *saved_mem, *saved_mem1;
            for (int i = 0; i < len; ++i) {
                query[i] = i < len / 2 ? 0 : 3;
                if (round == 1) query[i] = 0;
                if (round == 2 && i % 47 == 0) query[i] = 4;
            }
            mem_collect_intv(&opt, &bwt, len, query, reference);
            mem_collect_intv(&opt, &bwt, len, query, aux);
            equal_intervals(&reference->mem, &aux->mem);
            equal_intervals(&reference->mem1, &aux->mem1);
            nonempty += aux->mem.n > 0;
            assert((aux->tmpv[0] == &aux->ldm_tmpv[0]) ==
                   (aux->ldm_tmp_storage && len < SWBWA_SMEM_LDM_TMP_CAPACITY));
            if (aux->tmpv[0] == &aux->ldm_tmpv[0]) {
                assert(aux->ldm_tmpv[0].a == aux->ldm_tmp_storage);
                assert(aux->ldm_tmpv[1].a == aux->ldm_tmp_storage +
                       SWBWA_SMEM_LDM_TMP_CAPACITY);
                assert(aux->ldm_tmpv[0].m == SWBWA_SMEM_LDM_TMP_CAPACITY);
                assert(aux->ldm_tmpv[1].m == SWBWA_SMEM_LDM_TMP_CAPACITY);
            }
            saved_mem = malloc((aux->mem.n + 1) * sizeof(bwtintv_t));
            saved_mem1 = malloc((aux->mem1.n + 1) * sizeof(bwtintv_t));
            assert(saved_mem && saved_mem1);
            if (aux->mem.n)
                memcpy(saved_mem, aux->mem.a, aux->mem.n * sizeof(bwtintv_t));
            if (aux->mem1.n)
                memcpy(saved_mem1, aux->mem1.a, aux->mem1.n * sizeof(bwtintv_t));
            /* Overwrite the entire arena and grow seeds beyond its capacity.
             * Collection's exported values must survive until the next read. */
            test_boundaries(aux);
            test_growth(aux);
            if (aux->mem.n)
                assert(!memcmp(saved_mem, aux->mem.a,
                               aux->mem.n * sizeof(bwtintv_t)));
            if (aux->mem1.n)
                assert(!memcmp(saved_mem1, aux->mem1.a,
                               aux->mem1.n * sizeof(bwtintv_t)));
            free(saved_mem);
            free(saved_mem1);
            equal_intervals(&reference->mem, &aux->mem);
            equal_intervals(&reference->mem1, &aux->mem1);
        }
    }
    assert(nonempty > 0);
    free(bwt.bwt);
    smem_aux_destroy(reference);
}

static void test_context(int refuse)
{
    void *context;
    smem_aux_t *aux;
    long context_bytes;
    size_t before;
    refuse_ldm = refuse;
    swbwa_ldm_begin_batch();
    context = swbwa_ldm_alloc(sizeof(worker12_context_t), 4);
    context_bytes = context ? sizeof(worker12_context_t) : 0;
    aux = smem_aux_init(1);
    check_ownership(aux, context_bytes);
    if (is_overlay(aux)) assert(swbwa_ldm_refusals() == 0);
    printf("pool=%d reuse=%d capacity=%d refuse=%d overlay=%d held=%ld "
           "intv_bytes=%zu smem_bytes=%zu arena_bytes=%zu\n",
           SWBWA_TEST_POOL_ALLOC, SWBWA_CHAIN_REUSE_SMEM_LDM,
           SWBWA_SMEM_LDM_TMP_CAPACITY,
           refuse, is_overlay(aux), swbwa_ldm_outstanding(), sizeof(bwtintv_t),
           smem_bytes(), aux->chain_arena ? smem_chain_arena_capacity(aux) : 0);
    test_collection_phases(aux);
    check_ownership(aux, context_bytes);
    before = allocations;
    for (int i = 0; i < 1000; ++i) {
        void *seed;
        aux->chain_arena_used = 0;
        seed = smem_chain_alloc(aux, 4 * sizeof(mem_seed_t), 1);
        smem_chain_free(aux, seed);
    }
    assert(allocations - before == (aux->chain_arena ? 0 : 1000));
    smem_aux_destroy(aux);
    assert(swbwa_ldm_outstanding() == context_bytes);
    if (context) swbwa_ldm_release(context, context_bytes);
    assert(swbwa_ldm_outstanding() == 0);
    assert(swbwa_ldm_peak() <= (40 << 10));
}

int main(void)
{
    smem_aux_t *aux;
    void *seed;
#if SWBWA_TEST_POOL_ALLOC
    void *pool = malloc(16 << 20);
    assert(pool);
    set_big_buffer(pool, 16 << 20);
#endif
    test_context(0);
    test_context(1);
    /* Temporary aux must not lend storage to returned seeds. */
    refuse_ldm = 0;
    aux = smem_aux_init(0);
    assert(!aux->ldm_tmp_storage && !aux->chain_arena);
    seed = smem_chain_alloc(aux, sizeof(mem_seed_t), 1);
    smem_aux_destroy(aux);
    memset(seed, 0x5a, sizeof(mem_seed_t));
    smem_chain_free(NULL, seed);
#if SWBWA_TEST_POOL_ALLOC
    free(pool);
#endif
    puts("SMEM overlay collection/lifecycle checks passed");
    return 0;
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix="swbwa-smem-overlay-") as directory:
        directory = Path(directory)
        harness = directory / "test.c"
        harness.write_text(HARNESS)
        (directory / "crts.h").write_text(
            "void athread_lock(volatile long *);\n"
            "void athread_unlock(volatile long *);\n"
            "void athread_ssync_array(void);\n"
            "void athread_ssync_node(void);\n")
        # Report actual sizeof: default 256 means 16 KiB with 32-byte intervals.
        variants = [(capacity, reuse) for capacity in (256, 512) for reuse in (0, 1)]
        variants += [(384, 1), (64, 1), (0, 1)]  # Exact/small/no backing.
        for pool in (1, 0):
            for capacity, reuse in variants:
                executable = directory / f"test-{pool}-{capacity}-{reuse}"
                command = [os.environ.get("CC", "cc"), "-std=gnu11", "-O1", "-g",
                           "-fsanitize=address,undefined", "-fno-sanitize-recover=all",
                           "-fno-omit-frame-pointer",
                           "-ffunction-sections", "-fdata-sections",
                           "-Werror=implicit-function-declaration", "-D__uncached=",
                           "-D_PEN=swbwa_test_cpe_id",
                           "-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=0",
                           "-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=0",
                           f"-DSWBWA_CHAIN_REUSE_SMEM_LDM={reuse}",
                           f"-DSWBWA_SMEM_LDM_TMP_CAPACITY={capacity}",
                           f"-DSWBWA_TEST_POOL_ALLOC={pool}",
                           f"-I{directory}", "-Itests/stubs", "-Isrc/slave", "-Iinclude",
                           "-I.", "-include", "include/swbwa_config.h", str(harness),
                           "src/slave/malloc_wrap.c",
                           "-Wl,--gc-sections", "-lm", "-o", str(executable)]
                subprocess.run(command, cwd=ROOT, check=True)
                subprocess.run([str(executable)], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
