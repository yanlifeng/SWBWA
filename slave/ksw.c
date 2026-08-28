/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#if defined __SSE2__
#include <emmintrin.h>
#elif defined __ARM_NEON
#include "neon_sse.h"
#else
#include "scalar_sse.h"
#endif
#include "ksw.h"
#include "swbwa_cpe_profile.h"

#if SWBWA_ENABLE_CPE_MALLOC_WRAPPER
#  include "malloc_wrap.h"
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#define SWBWA_MATESW_HOT __attribute__((optimize("O3")))
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#define SWBWA_MATESW_HOT
#endif

#include <slave.h>

#if defined __SSE2__ || defined __ARM_NEON
#define SWBWA_KSW_I16_SHIFT_LEFT_ONE(value) _mm_slli_si128((value), 2)
#else
#define SWBWA_KSW_I16_SHIFT_LEFT_ONE(value) \
	swbwa_i16_shift_left_lane(value)
#endif

const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };

struct _kswq_t {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
	size_t allocation_bytes;
	int in_ldm;
};

enum { SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES = 64 << 10 };

#if SWBWA_ENABLE_CPE_PROFILE
enum {
	SWBWA_MATESW_KSW_FORWARD = 0,
	SWBWA_MATESW_KSW_REVERSE = 1
};

static __thread swbwa_matesw_ksw_work_t last_matesw_ksw_work;
static __thread swbwa_matesw_ksw_work_t *active_matesw_ksw_work;
static __thread int active_matesw_ksw_phase;

void swbwa_matesw_ksw_work_take(swbwa_matesw_ksw_work_t *work)
{
	*work = last_matesw_ksw_work;
}

static void record_matesw_ksw_work(int rows, int slen,
		uint64_t lazy_steps)
{
	swbwa_matesw_ksw_work_t *work = active_matesw_ksw_work;

	if (work == NULL) return;
	if (active_matesw_ksw_phase == SWBWA_MATESW_KSW_FORWARD) {
		work->forward_rows = rows;
		work->forward_main_steps = (uint64_t)rows * slen;
		work->forward_lazy_steps = lazy_steps;
	} else {
		work->reverse_rows = rows;
		work->reverse_main_steps = (uint64_t)rows * slen;
		work->reverse_lazy_steps = lazy_steps;
	}
}

static void finish_matesw_ksw_work(swbwa_matesw_ksw_work_t *work)
{
	last_matesw_ksw_work = *work;
	active_matesw_ksw_work = NULL;
}

static __thread swbwa_matesw_ksw_work_t last_matesw_pair_work[2];

void swbwa_matesw_ksw_pair_work_take(swbwa_matesw_ksw_work_t work[2])
{
	work[0] = last_matesw_pair_work[0];
	work[1] = last_matesw_pair_work[1];
}
#endif

#if SWBWA_ENABLE_FLOAT16_VECTOR && SWBWA_ENABLE_PACKED_INT8
#error "KSW FP16 and packed-int8 backends are mutually exclusive"
#endif

#if SWBWA_ENABLE_FLOAT16_VECTOR
enum { SWBWA_KSW_U8_LANES = SWBWA_KSW_U8_LOGICAL_LANES };
#elif SWBWA_ENABLE_PACKED_INT8
enum { SWBWA_KSW_U8_LANES = 32 };
#else
enum { SWBWA_KSW_U8_LANES = 16 };
#endif

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */


//__thread_local_fix char kswq_fix[32 << 10];


static kswq_t *ksw_qinit_impl(int size, int qlen, const uint8_t *query,
		int m, const int8_t *mat, int prefer_ldm)
{
	kswq_t *q;
	int slen, a, tmp, p;
	size_t allocation_bytes;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	if (size == 1) p = SWBWA_KSW_U8_LANES;
	slen = (qlen + p - 1) / p; // segmented length

	allocation_bytes = sizeof(kswq_t) + 63 +
	                   sizeof(__m128i) * slen * (m + 4);
	if (prefer_ldm &&
	    allocation_bytes <= SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES)
		q = (kswq_t*)ldm_malloc(allocation_bytes);
	else
		q = (kswq_t*)malloc(allocation_bytes);
	assert(q != NULL);
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 63) >> 6 << 6); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	q->allocation_bytes = allocation_bytes;
	q->in_ldm = prefer_ldm &&
	            allocation_bytes <= SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
#if SWBWA_ENABLE_PACKED_INT8
		int16_t *t = (int16_t*)q->qp;
#elif SWBWA_ENABLE_FLOAT16_VECTOR
		_Float16 *t = (_Float16*)q->qp;
#else
		int *t = (int*)q->qp;
#endif

		for (a = 0; a < m; ++a) {
			int i, k;
			const int8_t *ma = mat + a * m;
#if SWBWA_ENABLE_FLOAT16_VECTOR
			for (i = 0; i < slen; ++i) {
				int lane;

				for (lane = 0; lane < 32; ++lane) {
					k = i + lane * slen;
					*t++ = lane < p
					     ? (k >= qlen ? 0 : ma[query[k]]) + q->shift
					     : 0;
				}
			}
#else
			int nlen = slen * p;

			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
#endif
		}
	} else {
	#if !defined __SSE2__ && !defined __ARM_NEON
		int32_t *t = (int32_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i;
			const int8_t *ma = mat + a * m;

			for (i = 0; i < slen; ++i) {
				int lane;

				for (lane = 0; lane < 16; ++lane) {
					int k = i + lane * slen;

					*t++ = lane < p
					     ? (k >= qlen ? 0 : ma[query[k]])
					     : 0;
				}
			}
		}
	#else
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	#endif
	}
	return q;
}

kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m,
		const int8_t *mat)
{
	return ksw_qinit_impl(size, qlen, query, m, mat, 0);
}

static void ksw_qdestroy(kswq_t *q)
{
	if (q == NULL) return;
	if (q->in_ldm)
		ldm_free(q, q->allocation_bytes);
	else
		free(q);
}

#if defined __ARM_NEON
// This macro implicitly uses each function's `zero` local variable
#define _mm_slli_si128(a, n) (vextq_u8(zero, (a), 16 - (n)))
#endif

SWBWA_MATESW_HOT
kswr_t ksw_u8(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
	kswr_t r;
#if SWBWA_ENABLE_CPE_PROFILE
	uint64_t profile_lazy_steps = 0;
#endif

#if defined __SSE2__
#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

// Given entries with arbitrary values, return whether they are all 0x00
#define allzero_16(xx) (_mm_movemask_epi8(_mm_cmpeq_epi8((xx), zero)) == 0xffff)

#elif defined __ARM_NEON
#define __max_16(ret, xx) (ret) = vmaxvq_u8((xx))
#define allzero_16(xx) (vmaxvq_u8((xx)) == 0)

#else
#define __max_16(ret, xx) (ret) = m128i_max_u8((xx))
#define allzero_16(xx) (m128i_allzero((xx)))
#endif

	const int p_new = SWBWA_KSW_U8_LANES;

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi8(_o_del + _e_del);
	e_del = _mm_set1_epi8(_e_del);
	oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
	e_ins = _mm_set1_epi8(_e_ins);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
#if SWBWA_ENABLE_CPE_PROFILE
		uint64_t profile_row_lazy_steps =
			(uint64_t)SWBWA_KSW_U8_LANES * slen;
#endif
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, e_ins);
			t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < p_new); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, oe_ins);
				f = _mm_subs_epu8(f, e_ins);
				if (UNLIKELY(allzero_16(_mm_subs_epu8(f, h)))) {
#if SWBWA_ENABLE_CPE_PROFILE
					profile_row_lazy_steps = (uint64_t)k * slen + j + 1;
#endif
					goto end_loop16;
				}
			}
		}
end_loop16:
#if SWBWA_ENABLE_CPE_PROFILE
		profile_lazy_steps += profile_row_lazy_steps;
#endif
		//int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
#if SWBWA_ENABLE_CPE_PROFILE
	record_matesw_ksw_work(i < tlen ? i + 1 : tlen, slen,
	                       profile_lazy_steps);
#endif
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high;
#if !SWBWA_ENABLE_FLOAT16_VECTOR
		int qlen = slen * p_new;
#endif
		//uint8_t *t = (uint8_t*)Hmax;
#if SWBWA_ENABLE_PACKED_INT8
		uint16_t *t = (uint16_t*)Hmax;
#elif !SWBWA_ENABLE_FLOAT16_VECTOR
		int *t = (int*)Hmax;
#endif
#if SWBWA_ENABLE_FLOAT16_VECTOR
		for (i = 0; i < slen; ++i) {
			int lane;
			_Float16 *v = (_Float16*)(Hmax + i);

			for (lane = 0; lane < p_new; ++lane) {
				tmp = i + lane * slen;
				if ((int)v[lane] > max) max = v[lane], r.qe = tmp;
				else if ((int)v[lane] == max && tmp < r.qe) r.qe = tmp;
			}
		}
#else
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / p_new + i % p_new * slen;
			else if ((int)*t == max && (tmp = i / p_new + i % p_new * slen) < r.qe) r.qe = tmp; 
#endif
		//printf("%d,%d\n", max, gmax);
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

SWBWA_MATESW_HOT
kswr_t ksw_i16(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;

#if defined __SSE2__
#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)

// Given entries all either 0x0000 or 0xffff, return whether they are all 0x0000
#define allzero_0f_8(xx) (!_mm_movemask_epi8((xx)))

#elif defined __ARM_NEON
#define __max_8(ret, xx) (ret) = vmaxvq_s16(vreinterpretq_s16_u8((xx)))
#define allzero_0f_8(xx) (vmaxvq_u16(vreinterpretq_u16_u8((xx))) == 0)

#else
#define __max_8(ret, xx) (ret) = m128i_max_s16((xx))
#define allzero_0f_8(xx) (swbwa_i16_allzero((xx)))
#endif

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi16(_o_del + _e_del);
	e_del = _mm_set1_epi16(_e_del);
	oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
	e_ins = _mm_set1_epi16(_e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = SWBWA_KSW_I16_SHIFT_LEFT_ONE(h);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, _mm_load_si128(S++));
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, e_del);
			t = _mm_subs_epu16(h, oe_del);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, e_ins);
			t = _mm_subs_epu16(h, oe_ins);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = SWBWA_KSW_I16_SHIFT_LEFT_ONE(f);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, oe_ins);
				f = _mm_subs_epu16(f, e_ins);
				if(UNLIKELY(allzero_0f_8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	r.score = gmax; r.te = te;
	{
		int max = -1, low, high;
	#if !defined __SSE2__ && !defined __ARM_NEON
		for (i = 0, r.qe = -1; i < slen; ++i) {
			int lane;
			int values[16];

			simd_store(Hmax[i].words, values);
			for (lane = 0; lane < 8; ++lane) {
				int query_pos = i + lane * slen;
				int score = values[lane];

				if (score > max)
					max = score, r.qe = query_pos;
				else if (score == max && query_pos < r.qe)
					r.qe = query_pos;
			}
		}
	#else
		int tmp, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp; 
	#endif
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

#if SWBWA_ENABLE_MATESW_DUAL_FORWARD && \
    SWBWA_KSW_U8_MODE == SWBWA_KSW_U8_INT32_16
enum {
	SWBWA_KSW_PAIR_QUERY_LENGTH = 150,
	SWBWA_KSW_PAIR_LANES = 16
};

typedef union {
	float16v32 values;
	intv16 words;
} swbwa_ksw_pair_vec_t;

typedef struct {
	int slen;
	uint8_t shift;
	uint8_t max;
	float16v32 *qp;
	float16v32 *H0;
	float16v32 *H1;
	float16v32 *E;
	float16v32 *Hmax;
	size_t allocation_bytes;
	int in_ldm;
} swbwa_ksw_pair_q_t;

static const intv16 swbwa_ksw_pair_low_mask = {
	-1, -1, -1, -1, -1, -1, -1, -1,
	 0,  0,  0,  0,  0,  0,  0,  0
};

static const intv16 swbwa_ksw_pair_shift_mask = {
	-1, -1, -1, -1, -1, -1, -1, -1,
	-65536, -1, -1, -1, -1, -1, -1, -1
};

static inline float16v32 ksw_pair_blend(float16v32 old_value,
		float16v32 new_value,
		int use_low, int use_high)
{
	swbwa_ksw_pair_vec_t old_bits, new_bits, result;
	intv16 mask;
	intv16 all_bits = -1;

	if (use_low && use_high) return new_value;
	if (!use_low && !use_high) return old_value;
	mask = use_low ? swbwa_ksw_pair_low_mask
	               : simd_vxorw(swbwa_ksw_pair_low_mask, all_bits);
	old_bits.values = old_value;
	new_bits.values = new_value;
	result.words = simd_vbisw(simd_vandw(new_bits.words, mask),
	                          simd_vandw(old_bits.words,
	                                      simd_vxorw(mask, all_bits)));
	return result.values;
}

static inline float16v32 ksw_pair_shift(float16v32 value)
{
	swbwa_ksw_pair_vec_t shifted;

	shifted.values = value;
	shifted.words = simd_vandw(simd_sllx(shifted.words, 16),
	                           swbwa_ksw_pair_shift_mask);
	return shifted.values;
}

static inline int ksw_pair_zero_mask(float16v32 value)
{
	swbwa_ksw_pair_vec_t bits, low, high;
	intv16 all_bits = -1;

	bits.values = value;
	low.words = simd_vandw(bits.words, swbwa_ksw_pair_low_mask);
	high.words = simd_vandw(
		bits.words, simd_vxorw(swbwa_ksw_pair_low_mask, all_bits));
	return (simd_reduc_smaxh(low.values) == (_Float16)0 ? 1 : 0) |
	       (simd_reduc_smaxh(high.values) == (_Float16)0 ? 2 : 0);
}

static inline void ksw_pair_max_values(float16v32 value, int maxima[2])
{
	swbwa_ksw_pair_vec_t bits, low, high;
	intv16 all_bits = -1;

	bits.values = value;
	low.words = simd_vandw(bits.words, swbwa_ksw_pair_low_mask);
	high.words = simd_vandw(
		bits.words, simd_vxorw(swbwa_ksw_pair_low_mask, all_bits));
	maxima[0] = (int)simd_reduc_smaxh(low.values);
	maxima[1] = (int)simd_reduc_smaxh(high.values);
}

static swbwa_ksw_pair_q_t *ksw_qinit_u8_pair(int qlen,
		const uint8_t *query0,
		const uint8_t *query1, int m, const int8_t *mat, int prefer_ldm)
{
	swbwa_ksw_pair_q_t *q;
	size_t allocation_bytes;
	int tmp;
	int a, i, lane;
	_Float16 profile[32] __attribute__((aligned(64)));

	q = NULL;
	allocation_bytes = sizeof(*q) + 63 +
	                   sizeof(float16v32) *
	                   ((qlen + SWBWA_KSW_PAIR_LANES - 1) /
	                    SWBWA_KSW_PAIR_LANES) * (m + 4);
	if (prefer_ldm &&
	    allocation_bytes <= SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES)
		q = (swbwa_ksw_pair_q_t *)ldm_malloc(allocation_bytes);
	else
		q = (swbwa_ksw_pair_q_t *)malloc(allocation_bytes);
	assert(q != NULL);
	q->slen = (qlen + SWBWA_KSW_PAIR_LANES - 1) /
	          SWBWA_KSW_PAIR_LANES;
	q->allocation_bytes = allocation_bytes;
	q->in_ldm = prefer_ldm &&
	            allocation_bytes <= SWBWA_KSW_LDM_QUERY_PROFILE_MAX_BYTES;
	q->qp = (float16v32 *)(((uintptr_t)q + sizeof(*q) + 63) >> 6 << 6);
	q->H0 = q->qp + q->slen * m;
	q->H1 = q->H0 + q->slen;
	q->E = q->H1 + q->slen;
	q->Hmax = q->E + q->slen;
	tmp = m * m;
	q->shift = 127;
	q->max = 0;
	for (a = 0; a < tmp; ++a) {
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->max) q->max = mat[a];
	}
	q->shift = 256 - q->shift;

	for (a = 0; a < m; ++a) {
		const int8_t *ma = mat + a * m;

		for (i = 0; i < q->slen; ++i) {
			for (lane = 0; lane < SWBWA_KSW_PAIR_LANES; ++lane) {
				int k = i + lane * q->slen;

				profile[lane] =
					(_Float16)((k >= qlen ? 0 : ma[query0[k]]) + q->shift);
				profile[SWBWA_KSW_PAIR_LANES + lane] =
					(_Float16)((k >= qlen ? 0 : ma[query1[k]]) + q->shift);
			}
			simd_load(q->qp[a * q->slen + i], profile);
		}
	}
	return q;
}

static void ksw_pair_qdestroy(swbwa_ksw_pair_q_t *q)
{
	if (q->in_ldm)
		ldm_free(q, q->allocation_bytes);
	else
		free(q);
}

static void ksw_pair_append_peak(uint64_t **peaks, int *count, int *capacity,
		int row, int score)
{
	if (*count == 0 || (int32_t)(*peaks)[*count - 1] + 1 != row) {
		if (*count == *capacity) {
			*capacity = *capacity ? *capacity << 1 : 8;
			*peaks = (uint64_t *)realloc(*peaks,
			                              sizeof(**peaks) * *capacity);
		}
		(*peaks)[(*count)++] = (uint64_t)score << 32 | row;
	} else if ((int)((*peaks)[*count - 1] >> 32) < score) {
		(*peaks)[*count - 1] = (uint64_t)score << 32 | row;
	}
}

static void ksw_pair_finish_result(swbwa_ksw_pair_q_t *q, int which,
		int gmax, int te,
		uint64_t *peaks, int peak_count, kswr_t *result)
{
	int i;
	int best = -1;
	_Float16 values[32] __attribute__((aligned(64)));

	*result = g_defr;
	result->score = gmax + q->shift < 255 ? gmax : 255;
	result->te = te;
	if (result->score == 255) return;

	/* Match ksw_u8's physical-lane scan and tie order exactly. */
	for (i = 0; i < q->slen * SWBWA_KSW_U8_LANES; ++i) {
		int lane = i % SWBWA_KSW_U8_LANES;
		int query_pos = i / SWBWA_KSW_U8_LANES + lane * q->slen;
		int score;

		if (lane == 0)
			simd_store(q->Hmax[i / SWBWA_KSW_U8_LANES], values);
		score = (int)values[which * SWBWA_KSW_PAIR_LANES + lane];
		if (score > best || (score == best && query_pos < result->qe))
			best = score, result->qe = query_pos;
	}
	if (result->score != 255 && peak_count > 0) {
		int radius = (result->score + q->max - 1) / q->max;
		int range_low = te - radius;
		int range_high = te + radius;

		for (i = 0; i < peak_count; ++i) {
			int end = (int32_t)peaks[i];
			int score = peaks[i] >> 32;

			if ((end < range_low || end > range_high) &&
			    score > result->score2)
				result->score2 = score, result->te2 = end;
		}
	}
}

static void ksw_u8_dual_forward(swbwa_ksw_pair_q_t *q, int tlen,
		const uint8_t *target0, const uint8_t *target1, int o_del,
		int e_del, int o_ins, int e_ins, int xtra, kswr_t results[2]
#if SWBWA_ENABLE_CPE_PROFILE
		, swbwa_matesw_ksw_work_t work[2]
#endif
		)
{
	int i, j, k;
	int active[2] = { 1, 1 };
	int gmax[2] = { 0, 0 };
	int te[2] = { -1, -1 };
	int minsc = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
	int endsc = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
	int peak_count[2] = { 0, 0 };
	int peak_capacity[2] = { 0, 0 };
	uint64_t *peaks[2] = { NULL, NULL };
	float16v32 v_oe_del = (_Float16)(o_del + e_del);
	float16v32 v_e_del = (_Float16)e_del;
	float16v32 v_oe_ins = (_Float16)(o_ins + e_ins);
	float16v32 v_e_ins = (_Float16)e_ins;
	float16v32 v_shift = (_Float16)q->shift;
	float16v32 zero = (_Float16)0;
	float16v32 limit = (_Float16)255;
	float16v32 *H0 = q->H0, *H1 = q->H1, *E = q->E;
	int slen = q->slen;

	for (i = 0; i < slen; ++i) {
		H0[i] = zero;
		E[i] = zero;
		q->Hmax[i] = zero;
	}
	for (i = 0; i < tlen; ++i) {
		int row_active[2] = { active[0], active[1] };
		int lazy_done[2] = { !active[0], !active[1] };
		int maxima[2];
		float16v32 h, e, t, f, row_max;
		float16v32 *score_profile1 = q->qp + target1[i] * slen;

#if SWBWA_ENABLE_CPE_PROFILE
		uint64_t row_lazy_steps[2] = {
			(uint64_t)SWBWA_KSW_U8_LANES * slen,
			(uint64_t)SWBWA_KSW_U8_LANES * slen
		};
		for (j = 0; j < 2; ++j)
			if (row_active[j]) {
				++work[j].forward_rows;
				work[j].forward_main_steps += slen;
			}
#endif
		h = ksw_pair_shift(H0[slen - 1]);
		f = zero;
		row_max = zero;
		for (j = 0; j < slen; ++j) {
			float16v32 score = ksw_pair_blend(
				q->qp[target0[i] * slen + j], score_profile1[j], 0, 1);

			h = simd_sminh(simd_vaddh(h, score), limit);
			h = simd_smaxh(simd_vsubh(h, v_shift), zero);
			e = E[j];
			h = simd_smaxh(h, e);
			h = simd_smaxh(h, f);
			row_max = simd_smaxh(row_max, h);
			H1[j] = h;
			e = simd_smaxh(simd_vsubh(e, v_e_del), zero);
			t = simd_smaxh(simd_vsubh(h, v_oe_del), zero);
			e = simd_smaxh(e, t);
			E[j] = e;
			f = simd_smaxh(simd_vsubh(f, v_e_ins), zero);
			t = simd_smaxh(simd_vsubh(h, v_oe_ins), zero);
			f = simd_smaxh(f, t);
			h = H0[j];
		}
		for (k = 0; k < SWBWA_KSW_U8_LANES; ++k) {
			f = ksw_pair_shift(f);
			for (j = 0; j < slen; ++j) {
				int zero_mask;
				float16v32 previous_h;
				float16v32 updated_h;

				h = H1[j];
				previous_h = h;
				updated_h = simd_smaxh(previous_h, f);
				h = ksw_pair_blend(previous_h, updated_h,
				                   !lazy_done[0], !lazy_done[1]);
				H1[j] = h;
				h = simd_smaxh(simd_vsubh(h, v_oe_ins), zero);
				f = simd_smaxh(simd_vsubh(f, v_e_ins), zero);
				zero_mask = ksw_pair_zero_mask(
					simd_smaxh(simd_vsubh(f, h), zero));
				if (!lazy_done[0] && (zero_mask & 1)) {
					lazy_done[0] = 1;
#if SWBWA_ENABLE_CPE_PROFILE
					row_lazy_steps[0] = (uint64_t)k * slen + j + 1;
#endif
				}
				if (!lazy_done[1] && (zero_mask & 2)) {
					lazy_done[1] = 1;
#if SWBWA_ENABLE_CPE_PROFILE
					row_lazy_steps[1] = (uint64_t)k * slen + j + 1;
#endif
				}
				if (lazy_done[0] && lazy_done[1]) goto dual_lazy_done;
			}
		}
dual_lazy_done:
#if SWBWA_ENABLE_CPE_PROFILE
		for (j = 0; j < 2; ++j)
			if (row_active[j])
				work[j].forward_lazy_steps += row_lazy_steps[j];
#endif
		ksw_pair_max_values(row_max, maxima);
		for (j = 0; j < 2; ++j) {
			if (!row_active[j]) continue;
			if (maxima[j] >= minsc)
				ksw_pair_append_peak(&peaks[j], &peak_count[j],
				                     &peak_capacity[j], i, maxima[j]);
			if (maxima[j] > gmax[j]) {
				int p;

				gmax[j] = maxima[j];
				te[j] = i;
				for (p = 0; p < slen; ++p)
					q->Hmax[p] = ksw_pair_blend(
						q->Hmax[p], H1[p], j == 0, j == 1);
				if (gmax[j] + q->shift >= 255 || gmax[j] >= endsc)
					active[j] = 0;
			}
		}
		if (!active[0] && !active[1]) break;
		{
			float16v32 *swap = H1;
			H1 = H0;
			H0 = swap;
		}
	}
	for (i = 0; i < 2; ++i) {
		ksw_pair_finish_result(q, i, gmax[i], te[i], peaks[i],
		                       peak_count[i], &results[i]);
		free(peaks[i]);
	}
}
#endif

static inline void revseq(int l, uint8_t *s)
{
	int i, t;
	for (i = 0; i < l>>1; ++i)
		t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
}

static kswr_t ksw_align2_find_start(kswr_t r, int size, uint8_t *query,
		int tlen, uint8_t *target, int m, const int8_t *mat, int o_del,
		int e_del, int o_ins, int e_ins, int xtra, int profile_matesw)
{
	kswq_t *q;
	kswr_t rr;
	kswr_t (*func)(kswq_t*, int, const uint8_t*, int, int, int, int, int);

	if ((xtra & KSW_XSTART) == 0 ||
	    ((xtra & KSW_XSUBO) && r.score < (xtra & 0xffff)))
		return r;

	revseq(r.qe + 1, query);
	revseq(r.te + 1, target);
	if (profile_matesw)
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_REVERSE);
	q = ksw_qinit_impl(size, r.qe + 1, query, m, mat, profile_matesw);
#if SWBWA_ENABLE_CPE_PROFILE
	if (profile_matesw && active_matesw_ksw_work != NULL) {
		active_matesw_ksw_work->reverse_qlen = r.qe + 1;
		active_matesw_ksw_work->reverse_tlen = tlen;
		active_matesw_ksw_work->reverse_slen = q->slen;
		active_matesw_ksw_phase = SWBWA_MATESW_KSW_REVERSE;
	}
#endif
	if (profile_matesw)
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_REVERSE);

	func = q->size == 2 ? ksw_i16 : ksw_u8;
	if (profile_matesw)
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_DP_REVERSE);
	rr = func(q, tlen, target, o_del, e_del, o_ins, e_ins,
	          KSW_XSTOP | r.score);
	if (profile_matesw)
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_DP_REVERSE);
	revseq(r.qe + 1, query);
	revseq(r.te + 1, target);
	ksw_qdestroy(q);
	if (r.score == rr.score)
		r.tb = r.te - rr.te, r.qb = r.qe - rr.qe;
	return r;
}

static kswr_t ksw_align2_impl(int qlen, uint8_t *query, int tlen,
		uint8_t *target, int m, const int8_t *mat, int o_del, int e_del,
		int o_ins, int e_ins, int xtra, kswq_t **qry, int profile_matesw)
{
	int size;
	kswq_t *q;
	kswr_t r;
	kswr_t (*func)(kswq_t*, int, const uint8_t*, int, int, int, int, int);
#if SWBWA_ENABLE_CPE_PROFILE
	swbwa_matesw_ksw_work_t work;

	if (profile_matesw) {
		memset(&work, 0, sizeof(work));
		work.qlen = qlen;
		work.tlen = tlen;
		active_matesw_ksw_work = &work;
		active_matesw_ksw_phase = SWBWA_MATESW_KSW_FORWARD;
	}
#endif
	if (profile_matesw)
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_FORWARD);
	q = (qry && *qry)? *qry : ksw_qinit_impl((xtra&KSW_XBYTE)? 1 : 2,
			qlen, query, m, mat, profile_matesw && qry == 0);
	if (profile_matesw)
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_FORWARD);

	if (qry && *qry == 0) *qry = q;
	func = q->size == 2? ksw_i16 : ksw_u8;
	size = q->size;
#if SWBWA_ENABLE_CPE_PROFILE
	if (profile_matesw) {
		work.slen = q->slen;
		work.score_size = size;
	}
#endif
	if (profile_matesw)
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_DP_FORWARD);
	r = func(q, tlen, target, o_del, e_del, o_ins, e_ins, xtra);
	if (profile_matesw)
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_DP_FORWARD);
	if (qry == 0) ksw_qdestroy(q);
	//if (qry == 0) ldm_free(q, 32 << 10);

	r = ksw_align2_find_start(r, size, query, tlen, target, m, mat,
	                          o_del, e_del, o_ins, e_ins, xtra,
	                          profile_matesw);

#if SWBWA_ENABLE_CPE_PROFILE
	if (profile_matesw) finish_matesw_ksw_work(&work);
#endif

	return r;
}

kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target,
		int m, const int8_t *mat, int o_del, int e_del, int o_ins,
		int e_ins, int xtra, kswq_t **qry)
{
	return ksw_align2_impl(qlen, query, tlen, target, m, mat, o_del,
			e_del, o_ins, e_ins, xtra, qry, 0);
}

kswr_t ksw_align2_matesw(int qlen, uint8_t *query, int tlen,
		uint8_t *target, int m, const int8_t *mat, int o_del, int e_del,
		int o_ins, int e_ins, int xtra, kswq_t **qry)
{
	return ksw_align2_impl(qlen, query, tlen, target, m, mat, o_del,
			e_del, o_ins, e_ins, xtra, qry, 1);
}

void ksw_align2_matesw_dual_forward(
		int qlen0, uint8_t *query0, int tlen0, uint8_t *target0,
		int qlen1, uint8_t *query1, int tlen1, uint8_t *target1,
		int m, const int8_t *mat, int o_del, int e_del, int o_ins,
		int e_ins, int xtra, kswr_t results[2])
{
#if SWBWA_ENABLE_MATESW_DUAL_FORWARD && \
    SWBWA_KSW_U8_MODE == SWBWA_KSW_U8_INT32_16
	swbwa_ksw_pair_q_t *q;
#if SWBWA_ENABLE_CPE_PROFILE
	swbwa_matesw_ksw_work_t work[2];
#endif

	if ((xtra & KSW_XBYTE) && qlen0 == SWBWA_KSW_PAIR_QUERY_LENGTH &&
	    qlen0 == qlen1 && tlen0 == tlen1) {
#if SWBWA_ENABLE_CPE_PROFILE
		int i;

		memset(work, 0, sizeof(work));
		for (i = 0; i < 2; ++i) {
			work[i].qlen = qlen0;
			work[i].tlen = tlen0;
			work[i].slen = (qlen0 + SWBWA_KSW_PAIR_LANES - 1) /
			               SWBWA_KSW_PAIR_LANES;
			work[i].score_size = 1;
		}
#endif
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_FORWARD);
		q = ksw_qinit_u8_pair(qlen0, query0, query1, m, mat, 1);
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_QUERY_INIT_FORWARD);
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_KSW_DP_FORWARD);
		ksw_u8_dual_forward(q, tlen0, target0, target1, o_del, e_del,
		                    o_ins, e_ins, xtra, results
#if SWBWA_ENABLE_CPE_PROFILE
		                    , work
#endif
		                    );
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_KSW_DP_FORWARD);
		ksw_pair_qdestroy(q);

#if SWBWA_ENABLE_CPE_PROFILE
		active_matesw_ksw_work = &work[0];
#endif
		results[0] = ksw_align2_find_start(results[0], 1, query0, tlen0,
		                                    target0, m, mat, o_del, e_del,
		                                    o_ins, e_ins, xtra, 1);
#if SWBWA_ENABLE_CPE_PROFILE
		active_matesw_ksw_work = &work[1];
#endif
		results[1] = ksw_align2_find_start(results[1], 1, query1, tlen1,
		                                    target1, m, mat, o_del, e_del,
		                                    o_ins, e_ins, xtra, 1);
#if SWBWA_ENABLE_CPE_PROFILE
		last_matesw_pair_work[0] = work[0];
		last_matesw_pair_work[1] = work[1];
		active_matesw_ksw_work = NULL;
#endif
		return;
	}
#endif

	results[0] = ksw_align2_matesw(qlen0, query0, tlen0, target0, m,
	                                mat, o_del, e_del, o_ins, e_ins,
	                                xtra, NULL);
#if SWBWA_ENABLE_CPE_PROFILE
	last_matesw_pair_work[0] = last_matesw_ksw_work;
#endif
	results[1] = ksw_align2_matesw(qlen1, query1, tlen1, target1, m,
	                                mat, o_del, e_del, o_ins, e_ins,
	                                xtra, NULL);
#if SWBWA_ENABLE_CPE_PROFILE
	last_matesw_pair_work[1] = last_matesw_ksw_work;
#endif
}

kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry)
{
	return ksw_align2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, xtra, qry);
}

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

int ksw_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	assert(h0 > 0);
	// allocate memory
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; LIKELY(i < tlen); ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = beg; LIKELY(j < end); ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
			e = e > t? e : t;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m == 0) break;
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
		beg = j;
		for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
		end = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off)
{
	return ksw_extend2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
}

/********************
 * Global alignment *
 ********************/

#define MINUS_INF -0x40000000

static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

int ksw_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar_, uint32_t **cigar_)
{
	
	eh_t *eh;
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, score, n_col;
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	if (n_cigar_) *n_cigar_ = 0;
	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? malloc((long)n_col * tlen) : 0;
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = 0; eh[0].e = MINUS_INF;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	// DP loop
	for (i = 0; LIKELY(i < tlen); ++i) { // target sequence is in the outer loop
		int32_t f = MINUS_INF, h1, beg, end, t;
		int8_t *q = &qp[target[i] * qlen];
		beg = i > w? i - w : 0;
		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
		h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				e  = e > t? e : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1; k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1; // (i,k) points to the last cell
		while (i >= 0 && k >= 0) {
			which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}
	free(eh); free(qp); free(z);
	return score;
}

int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	return ksw_global2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, n_cigar_, cigar_);
}

/*******************************************
 * Main function (not compiled by default) *
 *******************************************/

#ifdef _KSW_MAIN

#include <unistd.h>
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, err_gzread)

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int main(int argc, char *argv[])
{
	int c, sa = 1, sb = 3, i, j, k, forward_only = 0, max_rseq = 0;
	int8_t mat[25];
	int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
	uint8_t *rseq = 0;
	gzFile fpt, fpq;
	kseq_t *kst, *ksq;

	// parse command line
	while ((c = getopt(argc, argv, "a:b:q:r:ft:1")) >= 0) {
		switch (c) {
			case 'a': sa = atoi(optarg); break;
			case 'b': sb = atoi(optarg); break;
			case 'q': gapo = atoi(optarg); break;
			case 'r': gape = atoi(optarg); break;
			case 't': minsc = atoi(optarg); break;
			case 'f': forward_only = 1; break;
			case '1': xtra |= KSW_XBYTE; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ksw [-1] [-f] [-a%d] [-b%d] [-q%d] [-r%d] [-t%d] <target.fa> <query.fa>\n", sa, sb, gapo, gape, minsc);
		return 1;
	}
	if (minsc > 0xffff) minsc = 0xffff;
	xtra |= KSW_XSUBO | minsc;
	// initialize scoring matrix
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
	// open file
	fpt = xzopen(argv[optind],   "r"); kst = kseq_init(fpt);
	fpq = xzopen(argv[optind+1], "r"); ksq = kseq_init(fpq);
	// all-pair alignment
	while (kseq_read(ksq) > 0) {
		kswq_t *q[2] = {0, 0};
		kswr_t r;
		for (i = 0; i < (int)ksq->seq.l; ++i) ksq->seq.s[i] = seq_nt4_table[(int)ksq->seq.s[i]];
		if (!forward_only) { // reverse
			if ((int)ksq->seq.m > max_rseq) {
				max_rseq = ksq->seq.m;
				rseq = (uint8_t*)realloc(rseq, max_rseq);
			}
			for (i = 0, j = ksq->seq.l - 1; i < (int)ksq->seq.l; ++i, --j)
				rseq[j] = ksq->seq.s[i] == 4? 4 : 3 - ksq->seq.s[i];
		}
		gzrewind(fpt); kseq_rewind(kst);
		while (kseq_read(kst) > 0) {
			for (i = 0; i < (int)kst->seq.l; ++i) kst->seq.s[i] = seq_nt4_table[(int)kst->seq.s[i]];
			r = ksw_align(ksq->seq.l, (uint8_t*)ksq->seq.s, kst->seq.l, (uint8_t*)kst->seq.s, 5, mat, gapo, gape, xtra, &q[0]);
			if (r.score >= minsc)
				err_printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", kst->name.s, r.tb, r.te+1, ksq->name.s, r.qb, r.qe+1, r.score, r.score2, r.te2);
			if (rseq) {
				r = ksw_align(ksq->seq.l, rseq, kst->seq.l, (uint8_t*)kst->seq.s, 5, mat, gapo, gape, xtra, &q[1]);
				if (r.score >= minsc)
					err_printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", kst->name.s, r.tb, r.te+1, ksq->name.s, (int)ksq->seq.l - r.qb, (int)ksq->seq.l - 1 - r.qe, r.score, r.score2, r.te2);
			}
		}
		free(q[0]); free(q[1]);
	}
	free(rseq);
	kseq_destroy(kst); err_gzclose(fpt);
	kseq_destroy(ksq); err_gzclose(fpq);
	return 0;
}
#endif
