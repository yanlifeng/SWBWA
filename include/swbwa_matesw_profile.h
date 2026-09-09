#ifndef SWBWA_MATESW_PROFILE_H
#define SWBWA_MATESW_PROFILE_H

#include <stdint.h>

enum {
    SWBWA_MATESW_CANDIDATE_BINS = 5,
    SWBWA_MATESW_RATIO_BINS = 5,
    SWBWA_MATESW_SAM_CANDIDATE_BINS = 6,
    SWBWA_MATESW_SAM_PROFILE_CAPACITY = 256
};

typedef struct {
    uint32_t qlen;
    uint32_t tlen;
    uint32_t slen;
    uint32_t score_size;
    uint32_t forward_rows;
    uint32_t reverse_qlen;
    uint32_t reverse_tlen;
    uint32_t reverse_slen;
    uint32_t reverse_rows;
    uint64_t forward_main_steps;
    uint64_t forward_lazy_steps;
    uint64_t reverse_main_steps;
    uint64_t reverse_lazy_steps;
} swbwa_matesw_ksw_work_t;

typedef struct {
    uint64_t invocations;
    uint64_t candidate_bins[SWBWA_MATESW_CANDIDATE_BINS];
    uint64_t candidates;
    uint64_t pairs;
    uint64_t paired_candidates;
    uint64_t same_orientation_pairs;
    uint64_t serial_work;
    uint64_t paired_work;
    uint64_t ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_pe_calls;
    uint64_t sam_pe_candidate_bins[SWBWA_MATESW_SAM_CANDIDATE_BINS];
    uint64_t sam_pe_both_directions;
    uint64_t same_read_pairs;
    uint64_t same_read_paired_candidates;
    uint64_t adjacent_read_groups;
    uint64_t adjacent_read_pairs;
    uint64_t adjacent_read_paired_candidates;
    uint64_t sam_u8_candidates;
    uint64_t sam_profiled_pairs;
    uint64_t sam_u8_pairs;
    uint64_t sam_qlen_equal_pairs;
    uint64_t sam_tlen_equal_pairs;
    uint64_t sam_forward_dimension_equal_pairs;
    uint64_t sam_reverse_dimension_equal_pairs;
    uint64_t sam_dimension_equal_pairs;
    uint64_t sam_profile_overflow;
    uint64_t sam_qlen_ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_tlen_ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_forward_work_ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_reverse_work_ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_total_work_ratio_bins[SWBWA_MATESW_RATIO_BINS];
    uint64_t sam_total_forward_main_steps;
    uint64_t sam_total_forward_lazy_steps;
    uint64_t sam_total_reverse_main_steps;
    uint64_t sam_total_reverse_lazy_steps;
    uint64_t sam_paired_forward_serial_work;
    uint64_t sam_paired_forward_lockstep_work;
    uint64_t sam_paired_reverse_serial_work;
    uint64_t sam_paired_reverse_lockstep_work;
} swbwa_matesw_profile_t;

#endif /* SWBWA_MATESW_PROFILE_H */
