#ifndef SCALAR_SSE_H
#define SCALAR_SSE_H

#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "simd.h"

#ifndef SWBWA_KSW_FUSED_GAP_UPDATE
#define SWBWA_KSW_FUSED_GAP_UPDATE SWBWA_ENABLE_CPE_KERNEL_OPT
#endif

#if SWBWA_KSW_FUSED_GAP_UPDATE != 0 && SWBWA_KSW_FUSED_GAP_UPDATE != 1
#error "SWBWA_KSW_FUSED_GAP_UPDATE must be 0 or 1"
#endif

#ifndef SWBWA_KSW_XOR_SELECT
#define SWBWA_KSW_XOR_SELECT SWBWA_ENABLE_CPE_KERNEL_OPT
#endif

#if SWBWA_KSW_XOR_SELECT != 0 && SWBWA_KSW_XOR_SELECT != 1
#error "SWBWA_KSW_XOR_SELECT must be 0 or 1"
#endif

#if SWBWA_KSW_XOR_SELECT
static inline intv16 swbwa_ksw_xor_select_words(intv16 a, intv16 b,
		intv16 choose_b)
{
	intv16 zero = 0;
	intv16 mask = simd_vsubw(zero, choose_b);

	/* Compare lanes are 0/1; expand before selecting all bits of a word. */
	return simd_vxorw(a, simd_vandw(simd_vxorw(a, b), mask));
}
#endif

# if SWBWA_ENABLE_FLOAT16_VECTOR

typedef union m128i {
    float16v32 val;
    intv16 words;
} __m128i;

#if SWBWA_KSW_U8_MODE == SWBWA_KSW_U8_FLOAT16_16
/* Keep the original KSW stripe width while using native FP16 arithmetic. */
static const intv16 swbwa_ksw_low_half_mask = {
    -1, -1, -1, -1, -1, -1, -1, -1,
     0,  0,  0,  0,  0,  0,  0,  0
};
#endif

static inline __m128i _mm_set1_epi32(int32_t n) {
	assert(n >= 0 && n <= 255);
	__m128i r;
    r.val = 1.0f * n;
    return r;
}

static inline __m128i _mm_load_si128(const __m128i *ptr) {
    __m128i r;
    simd_load(r.val, (_Float16 *)ptr);
    return r; 
}
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { 
    simd_store(a.val, (_Float16 *)ptr);
}

static inline int m128i_allzero(__m128i a) {
	return simd_reduc_smaxh(a.val) == (_Float16)0;
}

static inline __m128i _mm_slli_si128(__m128i a, int n) {
	__m128i r;
#if SWBWA_KSW_U8_MODE == SWBWA_KSW_U8_FLOAT16_16
	r.words = simd_vandw(simd_sllx(a.words, n * 16),
	                    swbwa_ksw_low_half_mask);
#else
	r.words = simd_sllx(a.words, n * 16);
#endif
	return r;
}

static inline __m128i _mm_max_epu8(__m128i a, __m128i b) {
    a.val = simd_smaxh(a.val, b.val);
	return a;
}

static inline __m128i _mm_min_epu8(__m128i a, __m128i b) {
    a.val = simd_sminh(a.val, b.val);
    return a;
}

static inline uint8_t m128i_max_u8(__m128i a) {
	return (uint8_t)simd_reduc_smaxh(a.val);
}

static inline __m128i _mm_set1_epi8(int8_t n) {
	__m128i r;
    r.val = 1.0f * n;
    return r;
}

static inline __m128i _mm_adds_epu8(__m128i a, __m128i b) {
    static float16v32 con_255 = 255.0f;
    a.val = simd_vaddh(a.val, b.val);
    a.val = simd_sminh(a.val, con_255);
	return a;
}

static inline __m128i _mm_subs_epu8(__m128i a, __m128i b) {
    static float16v32 con_0 = 0.0f;
    a.val = simd_vsubh(a.val, b.val);
    a.val = simd_smaxh(a.val, con_0);
	return a;
}

# else
typedef union m128i {
    intv16 val;
    intv16 words;
} __m128i;

intv16 v_min = 0;
intv16 v_max = -1;

static inline __m128i _mm_set1_epi32(int32_t n) {
	assert(n >= 0 && n <= 255);
	__m128i r;
    r.val = n;
    return r;
}

static inline __m128i _mm_load_si128(const __m128i *ptr) {
    __m128i r;
    simd_load(r.val, (int*)ptr); 
    return r; 
}
static inline void _mm_store_si128(__m128i *ptr, __m128i a) { 
    simd_store(a.val, (int *)ptr);
}

static inline int m128i_allzero(__m128i a) {
    intv16 xx = a.val;
    xx = simd_vbisw(xx, simd_sllx(xx, 8 * 32));
    xx = simd_vbisw(xx, simd_sllx(xx, 4 * 32));
    xx = simd_vbisw(xx, simd_sllx(xx, 2 * 32));
    xx = simd_vbisw(xx, simd_sllx(xx, 1 * 32));
    return simd_vextw15(xx) == 0;

}

static inline __m128i _mm_slli_si128(__m128i a, int n) {
    __m128i r;
    r.val = simd_sllx(a.val, n * 32);
	return r;
}

static inline intv16 _mm_max_intv16(intv16 a, intv16 b) {
#if SWBWA_KSW_XOR_SELECT
    return swbwa_ksw_xor_select_words(a, b, simd_vcmpltw(a, b));
#else
    intv16 mask = simd_vcmpltw(a, b);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 b_selected = simd_vandw(extended_mask, b);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 a_selected = simd_vandw(extended_mask, a);
    intv16 r = simd_vaddw(a_selected, b_selected);
	return r;
#endif
}

static inline intv16 _mm_min_intv16(intv16 a, intv16 b) {
#if SWBWA_KSW_XOR_SELECT
    return swbwa_ksw_xor_select_words(b, a, simd_vcmpltw(a, b));
#else
    intv16 mask = simd_vcmpltw(a, b);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 a_selected = simd_vandw(extended_mask, a);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 b_selected = simd_vandw(extended_mask, b);
    intv16 r = simd_vaddw(a_selected, b_selected);
	return r;
#endif
}

static inline __m128i _mm_max_epu8(__m128i a, __m128i b) {
#if SWBWA_KSW_XOR_SELECT
    a.val = _mm_max_intv16(a.val, b.val);
    return a;
#else
    
    intv16 mask = simd_vcmpltw(a.val, b.val);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 b_selected = simd_vandw(extended_mask, b.val);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 a_selected = simd_vandw(extended_mask, a.val);
    __m128i r;
    r.val = simd_vaddw(a_selected, b_selected);
	return r;
#endif
}

static inline __m128i _mm_min_epu8(__m128i a, __m128i b) {
#if SWBWA_KSW_XOR_SELECT
    a.val = _mm_min_intv16(a.val, b.val);
    return a;
#else
    
    intv16 mask = simd_vcmpltw(a.val, b.val);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 a_selected = simd_vandw(extended_mask, a.val);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 b_selected = simd_vandw(extended_mask, b.val);
    __m128i r;
    r.val = simd_vaddw(a_selected, b_selected);
	return r;
#endif
}

static inline uint8_t m128i_max_u8(__m128i a) {
    int val[16];
    simd_store(a.val, &(val[0]));
	int max = 0;
	for (int i = 0; i < 16; i++)
		if (max < val[i]) max = val[i];
	return max;
}

static inline __m128i _mm_set1_epi8(int8_t n) {
	__m128i r;
    r.val = n;
    return r;
}

static inline __m128i _mm_adds_epu8(__m128i a, __m128i b) {
    static intv16 con_255 = 255;
    a.val = simd_vaddw(a.val, b.val);
#if SWBWA_KSW_XOR_SELECT
    a.val = _mm_min_intv16(a.val, con_255);
#else
    intv16 mask = simd_vcmpltw(a.val, con_255);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 a_selected = simd_vandw(extended_mask, a.val);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 b_selected = simd_vandw(extended_mask, con_255);
    a.val = simd_vaddw(a_selected, b_selected);

#endif
	return a;
}
static inline __m128i _mm_subs_epu8(__m128i a, __m128i b) {
    static intv16 con_0 = 0;
    a.val = simd_vsubw(a.val, b.val);
#if SWBWA_KSW_XOR_SELECT
    a.val = simd_vandw(a.val,
                       simd_vsubw(con_0, simd_vcmpltw(con_0, a.val)));
#else
    intv16 mask = simd_vcmpltw(a.val, con_0);
    intv16 extended_mask = simd_vsubw(v_min, mask);
    intv16 b_selected = simd_vandw(extended_mask, con_0);
    extended_mask = simd_vxorw(extended_mask, v_max);
    intv16 a_selected = simd_vandw(extended_mask, a.val);
    a.val = simd_vaddw(a_selected, b_selected);
	
#endif
	return a;
}
# endif

static const intv16 swbwa_i16_low8_mask = {
	-1, -1, -1, -1, -1, -1, -1, -1,
	 0,  0,  0,  0,  0,  0,  0,  0
};

static inline intv16 swbwa_i16_select(intv16 a, intv16 b, intv16 choose_b)
{
#if SWBWA_KSW_XOR_SELECT
	return swbwa_ksw_xor_select_words(a, b, choose_b);
#else
	intv16 zero = 0;
	intv16 all_bits = -1;
	intv16 mask = simd_vsubw(zero, choose_b);

	return simd_vbisw(simd_vandw(mask, b),
	                  simd_vandw(simd_vxorw(mask, all_bits), a));
#endif
}

static inline intv16 swbwa_i16_max_words(intv16 a, intv16 b)
{
	return swbwa_i16_select(a, b, simd_vcmpltw(a, b));
}

static inline intv16 swbwa_i16_min_words(intv16 a, intv16 b)
{
	return swbwa_i16_select(b, a, simd_vcmpltw(a, b));
}

#if SWBWA_KSW_FUSED_GAP_UPDATE
/* max(max(gap-ext, 0), max(h-open_ext, 0)) needs only one zero clamp.
 * Keep signed word differences until the final clamp; no lane is repacked.
 */
static inline __m128i swbwa_ksw_gap_update_words(__m128i gap,
		__m128i ext, __m128i h, __m128i open_ext)
{
	intv16 zero = 0;
	intv16 value = swbwa_i16_max_words(simd_vsubw(gap.words, ext.words),
	                                  simd_vsubw(h.words, open_ext.words));
	__m128i r;

	r.words = simd_vandw(value,
	                    simd_vsubw(zero, simd_vcmpltw(zero, value)));
	return r;
}

static inline __m128i swbwa_ksw_i16_gap_update(__m128i gap,
		__m128i ext, __m128i h, __m128i open_ext)
{
	__m128i r = swbwa_ksw_gap_update_words(gap, ext, h, open_ext);

	r.words = simd_vandw(r.words, swbwa_i16_low8_mask);
	return r;
}
#endif

static inline __m128i swbwa_i16_shift_left_lane(__m128i a)
{
	__m128i r;

	r.words = simd_vandw(simd_sllx(a.words, 32), swbwa_i16_low8_mask);
	return r;
}

static inline int swbwa_i16_allzero(__m128i a)
{
	intv16 words = a.words;

	words = simd_vbisw(words, simd_sllx(words, 8 * 32));
	words = simd_vbisw(words, simd_sllx(words, 4 * 32));
	words = simd_vbisw(words, simd_sllx(words, 2 * 32));
	words = simd_vbisw(words, simd_sllx(words, 1 * 32));
	return simd_vextw15(words) == 0;
}

#if SWBWA_KSW_I16_MODE == SWBWA_KSW_I16_SCALAR_8
static inline __m128i _mm_adds_epi16(__m128i a, __m128i b)
{
	int av[16], bv[16];
	int i;

	simd_store(a.words, av);
	simd_store(b.words, bv);
	for (i = 0; i < 8; ++i) {
		int64_t sum = (int64_t)av[i] + bv[i];

		if (sum > 32767) sum = 32767;
		if (sum < -32768) sum = -32768;
		av[i] = (int)sum;
	}
	for (; i < 16; ++i) av[i] = 0;
	simd_load(a.words, av);
	return a;
}

static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b)
{
	int av[16], bv[16];
	int i;

	simd_store(a.words, av);
	simd_store(b.words, bv);
	for (i = 0; i < 8; ++i) av[i] = av[i] > bv[i] ? -1 : 0;
	for (; i < 16; ++i) av[i] = 0;
	simd_load(a.words, av);
	return a;
}

static inline __m128i _mm_max_epi16(__m128i a, __m128i b)
{
	int av[16], bv[16];
	int i;

	simd_store(a.words, av);
	simd_store(b.words, bv);
	for (i = 0; i < 8; ++i)
		if (av[i] < bv[i]) av[i] = bv[i];
	for (; i < 16; ++i) av[i] = 0;
	simd_load(a.words, av);
	return a;
}

static inline __m128i _mm_set1_epi16(int16_t n)
{
	int values[16] = { 0 };
	__m128i r;
	int i;

	for (i = 0; i < 8; ++i) values[i] = n;
	simd_load(r.words, values);
	return r;
}

static inline int16_t m128i_max_s16(__m128i a)
{
	int values[16];
	int max = -32768;
	int i;

	simd_store(a.words, values);
	for (i = 0; i < 8; ++i)
		if (max < values[i]) max = values[i];
	return (int16_t)max;
}

static inline __m128i _mm_subs_epu16(__m128i a, __m128i b)
{
	int av[16], bv[16];
	int i;

	simd_store(a.words, av);
	simd_store(b.words, bv);
	for (i = 0; i < 8; ++i) {
		int value = av[i] - bv[i];

		av[i] = value > 0 ? value : 0;
	}
	for (; i < 16; ++i) av[i] = 0;
	simd_load(a.words, av);
	return a;
}
#else
static inline __m128i _mm_adds_epi16(__m128i a, __m128i b)
{
	intv16 lower = -32768;
	intv16 upper = 32767;

	a.words = simd_vaddw(a.words, b.words);
	a.words = swbwa_i16_max_words(a.words, lower);
	a.words = swbwa_i16_min_words(a.words, upper);
	a.words = simd_vandw(a.words, swbwa_i16_low8_mask);
	return a;
}

static inline __m128i _mm_cmpgt_epi16(__m128i a, __m128i b)
{
	intv16 zero = 0;

	a.words = simd_vsubw(zero, simd_vcmpltw(b.words, a.words));
	a.words = simd_vandw(a.words, swbwa_i16_low8_mask);
	return a;
}

static inline __m128i _mm_max_epi16(__m128i a, __m128i b)
{
	a.words = swbwa_i16_max_words(a.words, b.words);
	a.words = simd_vandw(a.words, swbwa_i16_low8_mask);
	return a;
}

static inline __m128i _mm_set1_epi16(int16_t n)
{
	__m128i r;

	r.words = simd_vandw((intv16)n, swbwa_i16_low8_mask);
	return r;
}

static inline int16_t m128i_max_s16(__m128i a)
{
	int values[16];
	int max = -32768;
	int i;

	simd_store(a.words, values);
	for (i = 0; i < 8; ++i)
		if (max < values[i]) max = values[i];
	return (int16_t)max;
}

static inline __m128i _mm_subs_epu16(__m128i a, __m128i b)
{
	intv16 zero = 0;

	a.words = simd_vsubw(a.words, b.words);
	a.words = swbwa_i16_max_words(a.words, zero);
	a.words = simd_vandw(a.words, swbwa_i16_low8_mask);
	return a;
}
#endif

#endif
