/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "kstring.h"
#include "bwamem.h"
#include "swbwa_cpe_profile.h"
#include "kvec.h"
#include "utils.h"
#include "ksw.h"
#include <slave.h>

#if SWBWA_ENABLE_CPE_PROFILE
static __thread swbwa_matesw_profile_t local_matesw_profile;
static __thread uint64_t previous_sam_candidates;
static __thread int has_previous_sam;
typedef struct {
	uint32_t qlen;
	uint32_t tlen;
	uint32_t slen;
	uint32_t score_size;
	uint32_t reverse_qlen;
	uint32_t reverse_tlen;
	uint32_t reverse_slen;
	uint64_t forward_work;
	uint64_t reverse_work;
} swbwa_matesw_sam_candidate_t;
static __thread swbwa_matesw_sam_candidate_t
	direction0_candidates[SWBWA_MATESW_SAM_PROFILE_CAPACITY];
static __thread uint64_t direction_candidate_counts[2];
static __thread int active_sam_direction = -1;

static int matesw_ratio_bin(uint64_t first, uint64_t second)
{
	uint64_t larger = first > second ? first : second;
	uint64_t smaller = first > second ? second : first;

	if (larger == 0) return 0;
	if (smaller == 0) return SWBWA_MATESW_RATIO_BINS - 1;
	if (larger * 100 <= smaller * 110) return 0;
	if (larger * 100 <= smaller * 125) return 1;
	if (larger * 100 <= smaller * 150) return 2;
	if (larger <= smaller * 2) return 3;
	return 4;
}

void swbwa_matesw_profile_reset(void)
{
	memset(&local_matesw_profile, 0, sizeof(local_matesw_profile));
	previous_sam_candidates = 0;
	has_previous_sam = 0;
	direction_candidate_counts[0] = direction_candidate_counts[1] = 0;
	active_sam_direction = -1;
}

void swbwa_matesw_profile_commit(swbwa_matesw_profile_t *destination)
{
	int i;

	if (destination == NULL) return;
	destination->invocations += local_matesw_profile.invocations;
	destination->candidates += local_matesw_profile.candidates;
	destination->pairs += local_matesw_profile.pairs;
	destination->paired_candidates += local_matesw_profile.paired_candidates;
	destination->same_orientation_pairs +=
		local_matesw_profile.same_orientation_pairs;
	destination->serial_work += local_matesw_profile.serial_work;
	destination->paired_work += local_matesw_profile.paired_work;
	for (i = 0; i < SWBWA_MATESW_CANDIDATE_BINS; ++i)
		destination->candidate_bins[i] +=
			local_matesw_profile.candidate_bins[i];
	for (i = 0; i < SWBWA_MATESW_RATIO_BINS; ++i)
		destination->ratio_bins[i] += local_matesw_profile.ratio_bins[i];
	destination->sam_pe_calls += local_matesw_profile.sam_pe_calls;
	for (i = 0; i < SWBWA_MATESW_SAM_CANDIDATE_BINS; ++i)
		destination->sam_pe_candidate_bins[i] +=
			local_matesw_profile.sam_pe_candidate_bins[i];
	destination->sam_pe_both_directions +=
		local_matesw_profile.sam_pe_both_directions;
	destination->same_read_pairs += local_matesw_profile.same_read_pairs;
	destination->same_read_paired_candidates +=
		local_matesw_profile.same_read_paired_candidates;
	destination->adjacent_read_groups +=
		local_matesw_profile.adjacent_read_groups;
	destination->adjacent_read_pairs +=
		local_matesw_profile.adjacent_read_pairs;
	destination->adjacent_read_paired_candidates +=
		local_matesw_profile.adjacent_read_paired_candidates;
	destination->sam_u8_candidates += local_matesw_profile.sam_u8_candidates;
	destination->sam_profiled_pairs += local_matesw_profile.sam_profiled_pairs;
	destination->sam_u8_pairs += local_matesw_profile.sam_u8_pairs;
	destination->sam_qlen_equal_pairs +=
		local_matesw_profile.sam_qlen_equal_pairs;
	destination->sam_tlen_equal_pairs +=
		local_matesw_profile.sam_tlen_equal_pairs;
	destination->sam_forward_dimension_equal_pairs +=
		local_matesw_profile.sam_forward_dimension_equal_pairs;
	destination->sam_reverse_dimension_equal_pairs +=
		local_matesw_profile.sam_reverse_dimension_equal_pairs;
	destination->sam_dimension_equal_pairs +=
		local_matesw_profile.sam_dimension_equal_pairs;
	destination->sam_profile_overflow +=
		local_matesw_profile.sam_profile_overflow;
	for (i = 0; i < SWBWA_MATESW_RATIO_BINS; ++i) {
		destination->sam_qlen_ratio_bins[i] +=
			local_matesw_profile.sam_qlen_ratio_bins[i];
		destination->sam_tlen_ratio_bins[i] +=
			local_matesw_profile.sam_tlen_ratio_bins[i];
		destination->sam_forward_work_ratio_bins[i] +=
			local_matesw_profile.sam_forward_work_ratio_bins[i];
		destination->sam_reverse_work_ratio_bins[i] +=
			local_matesw_profile.sam_reverse_work_ratio_bins[i];
		destination->sam_total_work_ratio_bins[i] +=
			local_matesw_profile.sam_total_work_ratio_bins[i];
	}
	destination->sam_total_forward_main_steps +=
		local_matesw_profile.sam_total_forward_main_steps;
	destination->sam_total_forward_lazy_steps +=
		local_matesw_profile.sam_total_forward_lazy_steps;
	destination->sam_total_reverse_main_steps +=
		local_matesw_profile.sam_total_reverse_main_steps;
	destination->sam_total_reverse_lazy_steps +=
		local_matesw_profile.sam_total_reverse_lazy_steps;
	destination->sam_paired_forward_serial_work +=
		local_matesw_profile.sam_paired_forward_serial_work;
	destination->sam_paired_forward_lockstep_work +=
		local_matesw_profile.sam_paired_forward_lockstep_work;
	destination->sam_paired_reverse_serial_work +=
		local_matesw_profile.sam_paired_reverse_serial_work;
	destination->sam_paired_reverse_lockstep_work +=
		local_matesw_profile.sam_paired_reverse_lockstep_work;
}

static void begin_matesw_sam_profile(void)
{
	direction_candidate_counts[0] = direction_candidate_counts[1] = 0;
	active_sam_direction = -1;
}

static void set_matesw_sam_direction(int direction)
{
	assert(direction == 0 || direction == 1);
	active_sam_direction = direction;
}

static void record_matesw_sam_candidate(
	const swbwa_matesw_ksw_work_t *work)
{
	uint64_t forward_work = work->forward_main_steps +
	                        work->forward_lazy_steps;
	uint64_t reverse_work = work->reverse_main_steps +
	                        work->reverse_lazy_steps;
	uint64_t index;

	local_matesw_profile.sam_total_forward_main_steps +=
		work->forward_main_steps;
	local_matesw_profile.sam_total_forward_lazy_steps +=
		work->forward_lazy_steps;
	local_matesw_profile.sam_total_reverse_main_steps +=
		work->reverse_main_steps;
	local_matesw_profile.sam_total_reverse_lazy_steps +=
		work->reverse_lazy_steps;
	if (work->score_size == 1) ++local_matesw_profile.sam_u8_candidates;
	if (active_sam_direction < 0) return;

	index = direction_candidate_counts[active_sam_direction]++;
	if (active_sam_direction == 0) {
		swbwa_matesw_sam_candidate_t *candidate;

		if (index >= SWBWA_MATESW_SAM_PROFILE_CAPACITY) {
			++local_matesw_profile.sam_profile_overflow;
			return;
		}
		candidate = &direction0_candidates[index];
		candidate->qlen = work->qlen;
		candidate->tlen = work->tlen;
		candidate->slen = work->slen;
		candidate->score_size = work->score_size;
		candidate->reverse_qlen = work->reverse_qlen;
		candidate->reverse_tlen = work->reverse_tlen;
		candidate->reverse_slen = work->reverse_slen;
		candidate->forward_work = forward_work;
		candidate->reverse_work = reverse_work;
		return;
	}

	if (index < direction_candidate_counts[0]) {
		swbwa_matesw_sam_candidate_t *first;
		uint64_t first_total, second_total;

		if (index >= SWBWA_MATESW_SAM_PROFILE_CAPACITY) {
			++local_matesw_profile.sam_profile_overflow;
			return;
		}
		first = &direction0_candidates[index];
		first_total = first->forward_work + first->reverse_work;
		second_total = forward_work + reverse_work;
		++local_matesw_profile.sam_profiled_pairs;
		if (first->score_size == 1 && work->score_size == 1)
			++local_matesw_profile.sam_u8_pairs;
		if (first->qlen == work->qlen)
			++local_matesw_profile.sam_qlen_equal_pairs;
		if (first->tlen == work->tlen)
			++local_matesw_profile.sam_tlen_equal_pairs;
		if (first->qlen == work->qlen && first->tlen == work->tlen &&
		    first->slen == work->slen &&
		    first->score_size == work->score_size)
			++local_matesw_profile.sam_forward_dimension_equal_pairs;
		if (first->reverse_qlen == work->reverse_qlen &&
		    first->reverse_tlen == work->reverse_tlen &&
		    first->reverse_slen == work->reverse_slen)
			++local_matesw_profile.sam_reverse_dimension_equal_pairs;
		if (first->qlen == work->qlen && first->tlen == work->tlen &&
		    first->slen == work->slen &&
		    first->score_size == work->score_size &&
		    first->reverse_qlen == work->reverse_qlen &&
		    first->reverse_tlen == work->reverse_tlen &&
		    first->reverse_slen == work->reverse_slen)
			++local_matesw_profile.sam_dimension_equal_pairs;
		++local_matesw_profile.sam_qlen_ratio_bins[
			matesw_ratio_bin(first->qlen, work->qlen)];
		++local_matesw_profile.sam_tlen_ratio_bins[
			matesw_ratio_bin(first->tlen, work->tlen)];
		++local_matesw_profile.sam_forward_work_ratio_bins[
			matesw_ratio_bin(first->forward_work, forward_work)];
		++local_matesw_profile.sam_reverse_work_ratio_bins[
			matesw_ratio_bin(first->reverse_work, reverse_work)];
		++local_matesw_profile.sam_total_work_ratio_bins[
			matesw_ratio_bin(first_total, second_total)];
		local_matesw_profile.sam_paired_forward_serial_work +=
			first->forward_work + forward_work;
		local_matesw_profile.sam_paired_forward_lockstep_work +=
			first->forward_work > forward_work
			? first->forward_work : forward_work;
		local_matesw_profile.sam_paired_reverse_serial_work +=
			first->reverse_work + reverse_work;
		local_matesw_profile.sam_paired_reverse_lockstep_work +=
			first->reverse_work > reverse_work
			? first->reverse_work : reverse_work;
	}
}

static void end_matesw_sam_profile(void)
{
	active_sam_direction = -1;
}

static void record_matesw_pair_opportunity(
	int count, const int *work, const int *orientations)
{
	int order[4];
	int i;

	assert(count >= 0 && count <= 4);
	++local_matesw_profile.invocations;
	++local_matesw_profile.candidate_bins[count];
	local_matesw_profile.candidates += count;
	for (i = 0; i < count; ++i) {
		int j;

		order[i] = i;
		local_matesw_profile.serial_work += work[i];
		for (j = i; j > 0 && work[order[j - 1]] < work[order[j]]; --j) {
			int tmp = order[j - 1];
			order[j - 1] = order[j];
			order[j] = tmp;
		}
	}
	for (i = 0; i + 1 < count; i += 2) {
		uint64_t larger = work[order[i]];
		uint64_t smaller = work[order[i + 1]];

		++local_matesw_profile.pairs;
		local_matesw_profile.paired_candidates += 2;
		local_matesw_profile.paired_work += larger;
		if (orientations[order[i]] == orientations[order[i + 1]])
			++local_matesw_profile.same_orientation_pairs;
		if (larger * 100 <= smaller * 110)
			++local_matesw_profile.ratio_bins[0];
		else if (larger * 100 <= smaller * 125)
			++local_matesw_profile.ratio_bins[1];
		else if (larger * 100 <= smaller * 150)
			++local_matesw_profile.ratio_bins[2];
		else if (larger <= smaller * 2)
			++local_matesw_profile.ratio_bins[3];
		else
			++local_matesw_profile.ratio_bins[4];
	}
	if (i < count) local_matesw_profile.paired_work += work[order[i]];
}

static void record_matesw_sam_opportunity(const int candidates[2])
{
	uint64_t total = (uint64_t)candidates[0] + candidates[1];
	uint64_t pairs = candidates[0] < candidates[1]
	               ? candidates[0] : candidates[1];
	int bin = total < SWBWA_MATESW_SAM_CANDIDATE_BINS - 1
	        ? (int)total : SWBWA_MATESW_SAM_CANDIDATE_BINS - 1;

	++local_matesw_profile.sam_pe_calls;
	++local_matesw_profile.sam_pe_candidate_bins[bin];
	if (candidates[0] > 0 && candidates[1] > 0)
		++local_matesw_profile.sam_pe_both_directions;
	local_matesw_profile.same_read_pairs += pairs;
	local_matesw_profile.same_read_paired_candidates += pairs * 2;

	if (has_previous_sam) {
		pairs = previous_sam_candidates < total
		      ? previous_sam_candidates : total;
		++local_matesw_profile.adjacent_read_groups;
		local_matesw_profile.adjacent_read_pairs += pairs;
		local_matesw_profile.adjacent_read_paired_candidates += pairs * 2;
		has_previous_sam = 0;
	} else {
		previous_sam_candidates = total;
		has_previous_sam = 1;
	}
}
#else
static inline void record_matesw_pair_opportunity(
	int count, const int *work, const int *orientations)
{
	(void)count;
	(void)work;
	(void)orientations;
}

static inline void record_matesw_sam_opportunity(const int candidates[2])
{
	(void)candidates;
}
#endif

#if SWBWA_ENABLE_CPE_MALLOC_WRAPPER
#  include "malloc_wrap.h"
#endif


#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0
static int bwa_verbose = 1;

static inline int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist)
{
	int64_t p2;
	int r1 = (b1 >= l_pac), r2 = (b2 >= l_pac);
	p2 = r1 == r2? b2 : (l_pac<<1) - 1 - b2; // p2 is the coordinate of read 2 on the read 1 strand
	*dist = p2 > b1? p2 - b1 : b1 - p2;
	return (r1 == r2? 0 : 1) ^ (p2 > b1? 0 : 3);
}

static int cal_sub(const mem_opt_t *opt, mem_alnreg_v *r)
{
	int j;
	for (j = 1; j < r->n; ++j) { // choose unique alignment
		int b_max = r->a[j].qb > r->a[0].qb? r->a[j].qb : r->a[0].qb;
		int e_min = r->a[j].qe < r->a[0].qe? r->a[j].qe : r->a[0].qe;
		if (e_min > b_max) { // have overlap
			int min_l = r->a[j].qe - r->a[j].qb < r->a[0].qe - r->a[0].qb? r->a[j].qe - r->a[j].qb : r->a[0].qe - r->a[0].qb;
			if (e_min - b_max >= min_l * opt->mask_level) break; // significant overlap
		}
	}
	return j < r->n? r->a[j].score : opt->min_seed_len * opt->a;
}

__uncached long lock_f;

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int l_pos, int r_pos, const mem_alnreg_v *regs, mem_pestat_t pes[4], int* s_ids)
{
	int d, max;
	uint64_v isize[4];
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	memset(isize, 0, sizeof(kvec_t(int)) * 4);
	//for (int sid = 0; sid < n>>1; ++sid) {
	for (int sid = l_pos; sid < r_pos; ++sid) {
        //int i = s_ids[sid];
        int i = sid;
		int dir;
		int64_t is;
		mem_alnreg_v *r[2];
		r[0] = (mem_alnreg_v*)&regs[i<<1|0];
		r[1] = (mem_alnreg_v*)&regs[i<<1|1];
		if (r[0]->n == 0 || r[1]->n == 0) continue;
		if (cal_sub(opt, r[0]) > MIN_RATIO * r[0]->a[0].score) continue;
		if (cal_sub(opt, r[1]) > MIN_RATIO * r[1]->a[0].score) continue;
		if (r[0]->a[0].rid != r[1]->a[0].rid) continue; // not on the same chr
		dir = mem_infer_dir(l_pac, r[0]->a[0].rb, r[1]->a[0].rb, &is);
		if (is && is <= opt->max_ins) kv_push(uint64_t, isize[dir], is);
	}
	if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", __func__, isize[0].n, isize[1].n, isize[2].n, isize[3].n);
	for (d = 0; d < 4; ++d) { // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to merge these two.
		mem_pestat_t *r = &pes[d];
		uint64_v *q = &isize[d];
		int p25, p75, x, i;
		if (q->n < MIN_DIR_CNT) {
            //athread_lock(&lock_f);
			//printf("[M::%s] skip orientation %c%c as there are not enough pairs\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
            //athread_unlock(&lock_f);
			r->failed = 1;
			free(q->a);
			continue;
		} else {
            //athread_lock(&lock_f);
            //printf("[M::%s] analyzing insert size distribution for orientation %c%c...\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
            //athread_unlock(&lock_f);
        }
		ks_introsort_64(q->n, q->a);
		p25 = q->a[(int)(.25 * q->n + .499)];
		p75 = q->a[(int)(.75 * q->n + .499)];
		r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->low < 1) r->low = 1;
		r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
        //athread_lock(&lock_f);
		//printf("[M::%s] (25, 75) percentile: (%d, %d)\n", __func__, p25, p75);
		//printf("[M::%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, r->low, r->high);
        //athread_unlock(&lock_f);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
        //athread_lock(&lock_f);
		//printf("[M::%s] mean and std.dev: (%.2f, %.2f)\n", __func__, r->avg, r->std);
        //athread_unlock(&lock_f);
		r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->high < r->avg + MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->low < 1) r->low = 1;
        //athread_lock(&lock_f);
		//printf("[M::%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, r->low, r->high);
        //athread_unlock(&lock_f);
		free(q->a);
	}
	for (d = 0, max = 0; d < 4; ++d)
		max = max > isize[d].n? max : isize[d].n;
	for (d = 0; d < 4; ++d)
		if (pes[d].failed == 0 && isize[d].n < max * MIN_DIR_RATIO) {
			pes[d].failed = 1;
            //athread_lock(&lock_f);
			//printf("[M::%s] skip orientation %c%c\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
            //athread_unlock(&lock_f);
		}
}

typedef struct {
	uint8_t *seq;
	uint8_t *ref;
	int64_t rb;
	int64_t re;
	int is_rev;
} swbwa_matesw_candidate_t;

typedef struct {
	const mem_opt_t *opt;
	const bntseq_t *bns;
	const mem_alnreg_t *anchor;
	mem_alnreg_v *alignments;
	int l_ms;
	uint8_t *rev;
	int candidate_count;
	int added;
	int direction;
	swbwa_matesw_candidate_t candidates[4];
} swbwa_matesw_task_t;

static void swbwa_matesw_prepare(swbwa_matesw_task_t *task,
		const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
		const mem_pestat_t pes[4], const mem_alnreg_t *anchor, int l_ms,
		const uint8_t *mate_seq, mem_alnreg_v *alignments, int direction)
{
	int64_t l_pac = bns->l_pac;
	int i, r, skip[4];
#if SWBWA_ENABLE_CPE_PROFILE
	int candidate_work[4];
	int candidate_orientations[4];
#endif

	memset(task, 0, sizeof(*task));
	task->opt = opt;
	task->bns = bns;
	task->anchor = anchor;
	task->alignments = alignments;
	task->l_ms = l_ms;
	task->direction = direction;
	for (r = 0; r < 4; ++r) skip[r] = pes[r].failed ? 1 : 0;
	for (i = 0; i < alignments->n; ++i) {
		int64_t dist;

		r = mem_infer_dir(l_pac, anchor->rb, alignments->a[i].rb, &dist);
		if (dist >= pes[r].low && dist <= pes[r].high) skip[r] = 1;
	}
	for (r = 0; r < 4; ++r) {
		int is_rev, is_larger, rid = -1;
		uint8_t *seq, *ref = NULL;
		int64_t rb, re;
		swbwa_matesw_candidate_t *candidate;

		if (skip[r]) continue;
		is_rev = (r >> 1 != (r & 1));
		is_larger = !(r >> 1);
		if (is_rev) {
			if (task->rev == NULL) {
				task->rev = malloc(l_ms);
				for (i = 0; i < l_ms; ++i)
					task->rev[l_ms - 1 - i] =
						mate_seq[i] < 4 ? 3 - mate_seq[i] : 4;
			}
			seq = task->rev;
		} else {
			seq = (uint8_t *)mate_seq;
		}
		if (!is_rev) {
			rb = is_larger ? anchor->rb + pes[r].low
			               : anchor->rb - pes[r].high;
			re = (is_larger ? anchor->rb + pes[r].high
			                : anchor->rb - pes[r].low) + l_ms;
		} else {
			rb = (is_larger ? anchor->rb + pes[r].low
			                : anchor->rb - pes[r].high) - l_ms;
			re = is_larger ? anchor->rb + pes[r].high
			               : anchor->rb - pes[r].low;
		}
		if (rb < 0) rb = 0;
		if (re > l_pac << 1) re = l_pac << 1;
		if (rb < re) {
			swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_MATE_REF_FETCH);
			ref = bns_fetch_seq(bns, pac, &rb, (rb + re) >> 1, &re, &rid);
			swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_MATE_REF_FETCH);
		}
		if (anchor->rid != rid || re - rb < opt->min_seed_len) {
			free(ref);
			continue;
		}
		assert(task->candidate_count < 4);
		candidate = &task->candidates[task->candidate_count];
		candidate->seq = seq;
		candidate->ref = ref;
		candidate->rb = rb;
		candidate->re = re;
		candidate->is_rev = is_rev;
#if SWBWA_ENABLE_CPE_PROFILE
		candidate_work[task->candidate_count] =
			((l_ms + 15) >> 4) * (int)(re - rb);
		candidate_orientations[task->candidate_count] = is_rev;
#endif
		++task->candidate_count;
	}
#if SWBWA_ENABLE_CPE_PROFILE
	record_matesw_pair_opportunity(task->candidate_count, candidate_work,
	                                candidate_orientations);
#endif
}

static void swbwa_matesw_apply(swbwa_matesw_task_t *task, int index,
		kswr_t aln)
{
	swbwa_matesw_candidate_t *candidate = &task->candidates[index];
	const mem_alnreg_t *anchor = task->anchor;
	mem_alnreg_v *alignments = task->alignments;
	mem_alnreg_t b;
	int64_t l_pac = task->bns->l_pac;
	int i, insertion;

	if (aln.score < task->opt->min_seed_len || aln.qb < 0) return;
	memset(&b, 0, sizeof(b));
	b.rid = anchor->rid;
	b.is_alt = anchor->is_alt;
	b.qb = candidate->is_rev ? task->l_ms - (aln.qe + 1) : aln.qb;
	b.qe = candidate->is_rev ? task->l_ms - aln.qb : aln.qe + 1;
	b.rb = candidate->is_rev
	     ? (l_pac << 1) - (candidate->rb + aln.te + 1)
	     : candidate->rb + aln.tb;
	b.re = candidate->is_rev
	     ? (l_pac << 1) - (candidate->rb + aln.tb)
	     : candidate->rb + aln.te + 1;
	b.score = aln.score;
	b.csub = aln.score2;
	b.secondary = -1;
	b.seedcov = (b.re - b.rb < b.qe - b.qb ? b.re - b.rb : b.qe - b.qb) >> 1;
	kv_push(mem_alnreg_t, *alignments, b);
	for (i = 0; i < alignments->n - 1; ++i)
		if (alignments->a[i].score < b.score) break;
	insertion = i;
	for (i = alignments->n - 1; i > insertion; --i)
		alignments->a[i] = alignments->a[i - 1];
	alignments->a[insertion] = b;
	task->added = 1;
}

static void swbwa_matesw_run_one(swbwa_matesw_task_t *task, int index)
{
	swbwa_matesw_candidate_t *candidate = &task->candidates[index];
	const mem_opt_t *opt = task->opt;
	int xtra = KSW_XSUBO | KSW_XSTART |
	           (task->l_ms * opt->a < 250 ? KSW_XBYTE : 0) |
	           (opt->min_seed_len * opt->a);
	kswr_t aln;

#if SWBWA_ENABLE_CPE_PROFILE
	set_matesw_sam_direction(task->direction);
#endif
	swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_MATE_KSW_ALIGN);
	aln = ksw_align2_matesw(task->l_ms, candidate->seq,
	                        candidate->re - candidate->rb, candidate->ref,
	                        5, opt->mat, opt->o_del, opt->e_del, opt->o_ins,
	                        opt->e_ins, xtra, NULL);
	swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_MATE_KSW_ALIGN);
#if SWBWA_ENABLE_CPE_PROFILE
	{
		swbwa_matesw_ksw_work_t work;

		swbwa_matesw_ksw_work_take(&work);
		record_matesw_sam_candidate(&work);
	}
#endif
	swbwa_matesw_apply(task, index, aln);
}

static void swbwa_matesw_run_pair(swbwa_matesw_task_t tasks[2], int index)
{
	swbwa_matesw_candidate_t *first = &tasks[0].candidates[index];
	swbwa_matesw_candidate_t *second = &tasks[1].candidates[index];
	const mem_opt_t *opt = tasks[0].opt;
	int xtra = KSW_XSUBO | KSW_XSTART |
	           (tasks[0].l_ms * opt->a < 250 ? KSW_XBYTE : 0) |
	           (opt->min_seed_len * opt->a);
	kswr_t results[2];

	if (tasks[0].l_ms != tasks[1].l_ms) {
		swbwa_matesw_run_one(&tasks[0], index);
		swbwa_matesw_run_one(&tasks[1], index);
		return;
	}

	swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_MATE_KSW_ALIGN);
	ksw_align2_matesw_dual_forward(
		tasks[0].l_ms, first->seq, first->re - first->rb, first->ref,
		tasks[1].l_ms, second->seq, second->re - second->rb, second->ref,
		5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
		xtra, results);
	swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_MATE_KSW_ALIGN);
#if SWBWA_ENABLE_CPE_PROFILE
	{
		swbwa_matesw_ksw_work_t work[2];

		swbwa_matesw_ksw_pair_work_take(work);
		set_matesw_sam_direction(tasks[0].direction);
		record_matesw_sam_candidate(&work[0]);
		set_matesw_sam_direction(tasks[1].direction);
		record_matesw_sam_candidate(&work[1]);
	}
#endif
	swbwa_matesw_apply(&tasks[0], index, results[0]);
	swbwa_matesw_apply(&tasks[1], index, results[1]);
}

static void swbwa_matesw_finish(swbwa_matesw_task_t *task)
{
	extern int mem_sort_dedup_patch(const mem_opt_t *opt,
		const bntseq_t *bns, const uint8_t *pac, uint8_t *query,
		int n, mem_alnreg_t *a);
	int i;

	for (i = 0; i < task->candidate_count; ++i)
		free(task->candidates[i].ref);
	free(task->rev);
	if (task->added) {
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_MATE_DEDUP);
		task->alignments->n = mem_sort_dedup_patch(
			task->opt, NULL, NULL, NULL, task->alignments->n,
			task->alignments->a);
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_MATE_DEDUP);
	}
}

static int swbwa_mem_matesw_one(const mem_opt_t *opt, const bntseq_t *bns,
		const uint8_t *pac, const mem_pestat_t pes[4],
		const mem_alnreg_t *anchor, int l_ms, const uint8_t *mate_seq,
		mem_alnreg_v *alignments, int direction)
{
	swbwa_matesw_task_t task;
	int i;

	swbwa_matesw_prepare(&task, opt, bns, pac, pes, anchor, l_ms,
	                     mate_seq, alignments, direction);
	for (i = 0; i < task.candidate_count; ++i)
		swbwa_matesw_run_one(&task, i);
	swbwa_matesw_finish(&task);
	return task.candidate_count;
}

int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns,
		const uint8_t *pac, const mem_pestat_t pes[4],
		const mem_alnreg_t *anchor, int l_ms, const uint8_t *mate_seq,
		mem_alnreg_v *alignments)
{
	return swbwa_mem_matesw_one(opt, bns, pac, pes, anchor, l_ms,
	                            mate_seq, alignments, 0);
}

static int swbwa_mem_matesw_dual(const mem_opt_t *opt, const bntseq_t *bns,
		const uint8_t *pac, const mem_pestat_t pes[4],
		const mem_alnreg_t anchors[2], const int mate_lengths[2],
		const uint8_t *mate_seqs[2], mem_alnreg_v *alignments[2],
		int candidate_counts[2])
{
	swbwa_matesw_task_t tasks[2];
	int paired, i, total;

	for (i = 0; i < 2; ++i)
		swbwa_matesw_prepare(&tasks[i], opt, bns, pac, pes, &anchors[i],
		                     mate_lengths[i], mate_seqs[i], alignments[i], i);
	paired = tasks[0].candidate_count < tasks[1].candidate_count
	       ? tasks[0].candidate_count : tasks[1].candidate_count;
	for (i = 0; i < paired; ++i) swbwa_matesw_run_pair(tasks, i);
	for (i = paired; i < tasks[0].candidate_count; ++i)
		swbwa_matesw_run_one(&tasks[0], i);
	for (i = paired; i < tasks[1].candidate_count; ++i)
		swbwa_matesw_run_one(&tasks[1], i);
	total = tasks[0].candidate_count + tasks[1].candidate_count;
	if (candidate_counts != NULL) {
		candidate_counts[0] = tasks[0].candidate_count;
		candidate_counts[1] = tasks[1].candidate_count;
	}
	swbwa_matesw_finish(&tasks[0]);
	swbwa_matesw_finish(&tasks[1]);
	return total;
}

int mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v a[2], int id, int *sub, int *n_sub, int z[2], int n_pri[2])
{
	pair64_v v, u;
	int r, i, k, y[4], ret; // y[] keeps the last hit
	int64_t l_pac = bns->l_pac;
	kv_init(v); kv_init(u);
	for (r = 0; r < 2; ++r) { // loop through read number
		for (i = 0; i < n_pri[r]; ++i) {
			pair64_t key;
			mem_alnreg_t *e = &a[r].a[i];
			key.x = e->rb < l_pac? e->rb : (l_pac<<1) - 1 - e->rb; // forward position
			key.x = (uint64_t)e->rid<<32 | (key.x - bns->anns[e->rid].offset);
			key.y = (uint64_t)e->score << 32 | i << 2 | (e->rb >= l_pac)<<1 | r;
			kv_push(pair64_t, v, key);
		}
	}
	ks_introsort_128(v.n, v.a);
	y[0] = y[1] = y[2] = y[3] = -1;
	//for (i = 0; i < v.n; ++i) printf("[%d]\t%d\t%c%ld\n", i, (int)(v.a[i].y&1)+1, "+-"[v.a[i].y>>1&1], (long)v.a[i].x);
	for (i = 0; i < v.n; ++i) {
		for (r = 0; r < 2; ++r) { // loop through direction
			int dir = r<<1 | (v.a[i].y>>1&1), which;
			if (pes[dir].failed) continue; // invalid orientation
			which = r<<1 | ((v.a[i].y&1)^1);
			if (y[which] < 0) continue; // no previous hits
			for (k = y[which]; k >= 0; --k) { // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop takes a lot of time (I doubt)
				int64_t dist;
				int q;
				double ns;
				pair64_t *p;
				if ((v.a[k].y&3) != which) continue;
				dist = (int64_t)v.a[i].x - v.a[k].x;
				//printf("%d: %lld\n", k, dist);
				if (dist > pes[dir].high) break;
				if (dist < pes[dir].low)  continue;
				if(fabs(pes[dir].std) < 1e-6) {
                    q = 0;
                } else {
                    ns = (dist - pes[dir].avg) / pes[dir].std;
//                    double log_val, result;
//                    double erfc_val = erfc(fabs(ns) * M_SQRT1_2);
//                    double tmp = 2.0 * erfc_val;
//                    if (tmp <= 0.0) tmp = 1e-10;
//                    log_val = log(tmp);
//                    result = (v.a[i].y >> 32) + (v.a[k].y >> 32) + 0.721 * log_val * opt->a + 0.499;
//                    if (result > INT_MAX) result = INT_MAX;
//                    if (result < INT_MIN) result = INT_MIN;
//                    q = (int)llround(result);
                    q = (int)((v.a[i].y >> 32) + (v.a[k].y >> 32) + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * opt->a + 0.499); // .721 = 1/log(4)
                }
                if (q < 0) q = 0;
				p = kv_pushp(pair64_t, u);
				p->y = (uint64_t)k<<32 | i;
				p->x = (uint64_t)q<<32 | (hash_64(p->y ^ id<<8) & 0xffffffffU);
				//printf("[%lld,%lld]\t%d\tdist=%ld\n", v.a[k].x, v.a[i].x, q, (long)dist);
			}
		}
		y[v.a[i].y&3] = i;
	}
	if (u.n) { // found at least one proper pair
		int tmp = opt->a + opt->b;
		tmp = tmp > opt->o_del + opt->e_del? tmp : opt->o_del + opt->e_del;
		tmp = tmp > opt->o_ins + opt->e_ins? tmp : opt->o_ins + opt->e_ins;
		ks_introsort_128(u.n, u.a);
		i = u.a[u.n-1].y >> 32; k = u.a[u.n-1].y << 32 >> 32;
		z[v.a[i].y&1] = v.a[i].y<<32>>34; // index of the best pair
		z[v.a[k].y&1] = v.a[k].y<<32>>34;
		ret = u.a[u.n-1].x >> 32;
		*sub = u.n > 1? u.a[u.n-2].x>>32 : 0;
		for (i = (long)u.n - 2, *n_sub = 0; i >= 0; --i)
			if (*sub - (int)(u.a[i].x>>32) <= tmp) ++*n_sub;
	} else ret = 0, *sub = 0, *n_sub = 0;
	free(u.a); free(v.a);
	return ret;
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m);
void mem_reorder_primary5(int T, mem_alnreg_v *a);

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2])
{

	extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
	extern int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
	extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);
	extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const char *query);


	int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1, n_pri[2], n_aa[2];
	kstring_t str;
	mem_aln_t h[2], g[2], aa[2][2];
#if SWBWA_ENABLE_CPE_PROFILE
	int matesw_candidates[2] = { 0, 0 };
#endif

	str.l = str.m = 0; str.s = 0;
	memset(h, 0, sizeof(mem_aln_t) * 2);
	memset(g, 0, sizeof(mem_aln_t) * 2);
	n_aa[0] = n_aa[1] = 0;
	if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
		mem_alnreg_v b[2];
		swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_MATE_RESCUE);
#if SWBWA_ENABLE_CPE_PROFILE
		begin_matesw_sam_profile();
#endif
		kv_init(b[0]); kv_init(b[1]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < a[i].n; ++j)
				if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
					kv_push(mem_alnreg_t, b[i], a[i].a[j]);
#if SWBWA_ENABLE_MATESW_DUAL_FORWARD
		{
			int call_counts[2] = {
				b[0].n < opt->max_matesw ? b[0].n : opt->max_matesw,
				b[1].n < opt->max_matesw ? b[1].n : opt->max_matesw
			};
			int paired_calls = call_counts[0] < call_counts[1]
			                 ? call_counts[0] : call_counts[1];

			for (j = 0; j < paired_calls; ++j) {
				mem_alnreg_t anchors[2] = { b[0].a[j], b[1].a[j] };
				int mate_lengths[2] = { s[1].l_seq, s[0].l_seq };
				const uint8_t *mate_seqs[2] = {
					(uint8_t *)s[1].seq, (uint8_t *)s[0].seq
				};
				mem_alnreg_v *mate_alignments[2] = { &a[1], &a[0] };
				int candidate_counts[2];

				n += swbwa_mem_matesw_dual(
					opt, bns, pac, pes, anchors, mate_lengths, mate_seqs,
					mate_alignments, candidate_counts);
#if SWBWA_ENABLE_CPE_PROFILE
				matesw_candidates[0] += candidate_counts[0];
				matesw_candidates[1] += candidate_counts[1];
#endif
			}
			for (i = 0; i < 2; ++i) {
				for (j = paired_calls; j < call_counts[i]; ++j) {
					int count = swbwa_mem_matesw_one(
						opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq,
						(uint8_t *)s[!i].seq, &a[!i], i);

					n += count;
#if SWBWA_ENABLE_CPE_PROFILE
					matesw_candidates[i] += count;
#endif
				}
			}
		}
#else
		for (i = 0; i < 2; ++i) {
			for (j = 0; j < b[i].n && j < opt->max_matesw; ++j) {
				int count = swbwa_mem_matesw_one(
					opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq,
					(uint8_t *)s[!i].seq, &a[!i], i);

				n += count;
#if SWBWA_ENABLE_CPE_PROFILE
				matesw_candidates[i] += count;
#endif
			}
		}
#endif
		free(b[0].a); free(b[1].a);
#if SWBWA_ENABLE_CPE_PROFILE
		end_matesw_sam_profile();
#endif
		swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_MATE_RESCUE);
	}
#if SWBWA_ENABLE_CPE_PROFILE
	record_matesw_sam_opportunity(matesw_candidates);
#endif
	n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
	n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);
	if (opt->flag & MEM_F_PRIMARY5) {
		mem_reorder_primary5(opt->T, &a[0]);
		mem_reorder_primary5(opt->T, &a[1]);
	}


	if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;

	// pairing single-end hits
	//if (n_pri[0] && n_pri[1] && (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0) {
	if (n_pri[0] && n_pri[1]) {
        swbwa_cpe_profile_start(SWBWA_CPE_PROFILE_PAIRING);
        o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri);
        swbwa_cpe_profile_stop(SWBWA_CPE_PROFILE_PAIRING);
        if(o <= 0) goto no_pairing;
		int is_multi[2], q_pe, score_un, q_se[2];
		char **XA[2];
		// check if an end has multiple hits even after mate-SW
		for (i = 0; i < 2; ++i) {
			for (j = 1; j < n_pri[i]; ++j)
				if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
			is_multi[i] = j < n_pri[i]? 1 : 0;
		}
		if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
		// compute mapQ for the best SE hit
		score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
		//q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
		subo = subo > score_un? subo : score_un;
		q_pe = raw_mapq(o - subo, opt->a);
		if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
		if (q_pe < 0) q_pe = 0;
		if (q_pe > 60) q_pe = 60;
		q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);
		// the following assumes no split hits
		if (o > score_un) { // paired alignment is preferred
			mem_alnreg_t *c[2];
			c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
			for (i = 0; i < 2; ++i) {
				if (c[i]->secondary >= 0)
					c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
				q_se[i] = mem_approx_mapq_se(opt, c[i]);
			}
			q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
			q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
			extra_flag |= 2;
			// cap at the tandem repeat score
			q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
			q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
		} else { // the unpaired alignment is preferred
			z[0] = z[1] = 0;
			q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
			q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
		}
		for (i = 0; i < 2; ++i) {
			int k = a[i].a[z[i]].secondary_all;
			if (k >= 0 && k < n_pri[i]) { // switch secondary and primary if both of them are non-ALT
				assert(a[i].a[k].secondary_all < 0);
				for (j = 0; j < a[i].n; ++j)
					if (a[i].a[j].secondary_all == k || j == k)
						a[i].a[j].secondary_all = z[i];
				a[i].a[z[i]].secondary_all = -1;
			}
		}
		if (!(opt->flag & MEM_F_ALL)) {
			for (i = 0; i < 2; ++i)
				XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq);
		} else XA[0] = XA[1] = 0;
		// write SAM
		for (i = 0; i < 2; ++i) {
			h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]]);
			h[i].mapq = q_se[i];
			h[i].flag |= 0x40<<i | extra_flag;
			h[i].XA = XA[i]? XA[i][z[i]] : 0;
			aa[i][n_aa[i]++] = h[i];
			if (n_pri[i] < a[i].n) { // the read has ALT hits
				mem_alnreg_t *p = &a[i].a[n_pri[i]];
				if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
				g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
				g[i].flag |= 0x800 | 0x40<<i | extra_flag;
				g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
				aa[i][n_aa[i]++] = g[i];
			}
		}
		for (i = 0; i < n_aa[0]; ++i)
			mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]); // write read1 hits
		s[0].sam = strdup(str.s); str.l = 0;
		for (i = 0; i < n_aa[1]; ++i)
			mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]); // write read2 hits
		s[1].sam = str.s;
		if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
		// free
		for (i = 0; i < 2; ++i) {
			free(h[i].cigar); free(g[i].cigar);
			if (XA[i] == 0) continue;
			for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
			free(XA[i]);
		}
	} else goto no_pairing;

	return n;

no_pairing:

	for (i = 0; i < 2; ++i) {
		int which = -1;
		if (a[i].n) {
			if (a[i].a[0].score >= opt->T) which = 0;
			else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
				which = n_pri[i];
		}
		if (which >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which]);
		else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
	}
	if (!(opt->flag & MEM_F_NOPAIRING) && h[0].rid == h[1].rid && h[0].rid >= 0) { // if the top hits from the two ends constitute a proper pair, flag it.
		int64_t dist;
		int d;
		d = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
		if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
	}
	mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1]);
	mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0]);
	if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
	free(h[0].cigar); free(h[1].cigar);

	return n;
}
