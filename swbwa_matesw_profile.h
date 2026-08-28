#ifndef SWBWA_MATESW_PROFILE_H
#define SWBWA_MATESW_PROFILE_H

#include <stdint.h>

enum {
    SWBWA_MATESW_CANDIDATE_BINS = 5,
    SWBWA_MATESW_RATIO_BINS = 5
};

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
} swbwa_matesw_profile_t;

#endif /* SWBWA_MATESW_PROFILE_H */
