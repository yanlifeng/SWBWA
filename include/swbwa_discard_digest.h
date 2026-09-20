#ifndef SWBWA_DISCARD_DIGEST_H
#define SWBWA_DISCARD_DIGEST_H

#include <stddef.h>
#include <stdint.h>
#include "swbwa_config.h"

#ifndef SWBWA_CPE_DISCARD_DIGEST
#define SWBWA_CPE_DISCARD_DIGEST 0
#endif
#if SWBWA_CPE_DISCARD_DIGEST != 0 && SWBWA_CPE_DISCARD_DIGEST != 1
#error "SWBWA_CPE_DISCARD_DIGEST must be 0 or 1"
#endif
#define SWBWA_CPE_DISCARD_DIGEST_ACTIVE \
    (SWBWA_CPE_DISCARD_DIGEST && SWBWA_OUTPUT_MODE == SWBWA_OUTPUT_DISCARD)

#if SWBWA_CPE_DISCARD_DIGEST_ACTIVE
#if SWBWA_USE_MPI || !SWBWA_USE_CROSS_SEGMENT || \
    SWBWA_CPE_ALLOC_MODE != SWBWA_CPE_ALLOC_POOL
#error "CPE discard digest requires non-MPI cgs_cross+pool"
#endif
#if SWBWA_DISCARD_HASH_BYTES != 0
#error "CPE discard digest requires FULL hashing: SWBWA_DISCARD_HASH_BYTES=0"
#endif

#define SWBWA_DIGEST_RECORD_BYTES 24
#define SWBWA_DIGEST_NEEDS_HASH UINT64_C(1)
#define SWBWA_DIGEST_NEEDS_COUNT UINT64_C(2)
#define SWBWA_DIGEST_HASH_READY UINT64_C(3)
#define SWBWA_DIGEST_COUNT_READY UINT64_C(4)

static inline uint64_t swbwa_digest_load64(const unsigned char *p)
{
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

static inline void swbwa_digest_store64(unsigned char *p, uint64_t value)
{
    for (unsigned int i = 0; i < 8; ++i) p[i] = (unsigned char)(value >> (8 * i));
}

static inline uint64_t swbwa_digest_mix_word(const unsigned char *p)
{
    const uint64_t multiplier = UINT64_C(0xc6a4a7935bd1e995);
    uint64_t word = swbwa_digest_load64(p) * multiplier;
    word ^= word >> 47;
    return word * multiplier;
}

/* Exactly one existing output-write blob, excluding its terminating NUL. */
static inline uint64_t swbwa_digest_hash_blob(const void *data, size_t length)
{
    const uint64_t multiplier = UINT64_C(0xc6a4a7935bd1e995);
    const unsigned char *p = (const unsigned char *)data;
    size_t remaining = length;
    uint64_t hash = UINT64_C(0x9e3779b97f4a7c15) ^ ((uint64_t)length * multiplier);

    while (remaining >= 32) {
        uint64_t a = swbwa_digest_mix_word(p);
        uint64_t b = swbwa_digest_mix_word(p + 8);
        uint64_t c = swbwa_digest_mix_word(p + 16);
        uint64_t d = swbwa_digest_mix_word(p + 24);
        hash = (hash ^ a) * multiplier;
        hash = (hash ^ b) * multiplier;
        hash = (hash ^ c) * multiplier;
        hash = (hash ^ d) * multiplier;
        p += 32;
        remaining -= 32;
    }
    while (remaining >= 8) {
        hash = (hash ^ swbwa_digest_mix_word(p)) * multiplier;
        p += 8;
        remaining -= 8;
    }
    if (remaining) {
        uint64_t tail = 0;
        for (size_t i = 0; i < remaining; ++i) tail |= (uint64_t)p[i] << (8 * i);
        hash = (hash ^ tail) * multiplier;
    }
    hash ^= hash >> 47;
    hash *= multiplier;
    return hash ^ (hash >> 47);
}

/* The existing batch SAM ring owns these records through pipeline stage 3. */
static inline void swbwa_digest_prepare(void *record, uint64_t length, int enabled)
{
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && \
    __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    /* Preserve the byte-defined format; only aligned metadata uses word stores. */
    if (((uintptr_t)record & (uintptr_t)7) == 0) {
        unsigned char *aligned = (unsigned char *)__builtin_assume_aligned(record, 8);
        const uint64_t zero = 0;
        const uint64_t state = enabled ? SWBWA_DIGEST_NEEDS_HASH : SWBWA_DIGEST_NEEDS_COUNT;
        __builtin_memcpy(aligned, &length, sizeof(length));
        __builtin_memcpy(aligned + 8, &zero, sizeof(zero));
        __builtin_memcpy(aligned + 16, &state, sizeof(state));
        return;
    }
#endif
    unsigned char *p = (unsigned char *)record;
    swbwa_digest_store64(p, length);
    swbwa_digest_store64(p + 8, 0);
    swbwa_digest_store64(p + 16, enabled ? SWBWA_DIGEST_NEEDS_HASH : SWBWA_DIGEST_NEEDS_COUNT);
}

static inline void swbwa_digest_finish(void *record, const void *blob, size_t length)
{
    unsigned char *p = (unsigned char *)record;
    uint64_t state = swbwa_digest_load64(p + 16);
    uint64_t hash = 0;

    if (swbwa_digest_load64(p) != (uint64_t)length ||
        (state != SWBWA_DIGEST_NEEDS_HASH && state != SWBWA_DIGEST_NEEDS_COUNT) ||
        (length != 0 && blob == NULL)) {
        swbwa_digest_store64(p + 16, 0);
        return;
    }
    if (state == SWBWA_DIGEST_NEEDS_HASH && length != 0)
        hash = swbwa_digest_hash_blob(blob, length);
    swbwa_digest_store64(p + 8, hash);
    swbwa_digest_store64(p + 16, state == SWBWA_DIGEST_NEEDS_HASH ?
                          SWBWA_DIGEST_HASH_READY : SWBWA_DIGEST_COUNT_READY);
}

#ifdef __cplusplus
extern "C" {
#endif
int swbwa_output_discard_hash_active(void);
int swbwa_output_write_digest(const void *record);
#ifdef __cplusplus
}
#endif
#endif /* SWBWA_CPE_DISCARD_DIGEST_ACTIVE */
#endif /* SWBWA_DISCARD_DIGEST_H */
