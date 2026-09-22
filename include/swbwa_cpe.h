#ifndef SWBWA_CPE_H
#define SWBWA_CPE_H

#include "swbwa_config.h"
#include "swbwa_matesw_profile.h"
#include "swbwa_ldm_alloc.h"

#define SWBWA_CPE_CSR_COPY_BYTES (2UL << 20)
#define SWBWA_CPE_PRIVATE_BASE   0x400000000000UL

typedef struct {
    char *buffer;
    long long bytes_per_cpe;
} swbwa_cpe_pool_params_t;

/*
 * Shared MPE/CPE task parameters. Keep this structure pointer-only where
 * possible: it is passed directly between the two architectures.
 */
typedef struct {
    long work_item_count;
    void *worker_data;
    int *real_sizes;
    mem_alnreg_v *alignment_regions;
    char **sam_records;
    int *sam_lengths;
    const mem_pestat_t *pes;
    int *sequence_ids;
    long *profile_counters;
    swbwa_matesw_profile_t *matesw_profile;

    volatile int *completion_flags;
    volatile long *error_info; /* Cross runtime: CPE_COUNT slots, ERROR_WORDS each. */
    void *relocated_gp;
    unsigned long segment_offset;
    void *private_segment_copies;
    int tls_relocation_count;
    unsigned long *tls_relocations;

    char *fastq_buffer[2];
    char *formatted_buffer[2];
    long long fastq_size[2];
    long long formatted_buffer_size;
    long formatted_read_counts[SWBWA_CPE_COUNT];
    /* Bytes handed out by the CPE bump allocator, sampled per batch. */
    long pool_high_water[SWBWA_CPE_COUNT];
    /* LDM bytes still held after a batch (should be zero) and the batch peak. */
    long ldm_outstanding[SWBWA_CPE_COUNT];
    long ldm_peak[SWBWA_CPE_COUNT];
    long ldm_refusals[SWBWA_CPE_COUNT];
#if SWBWA_CPE_LDM_ALLOC
    swbwa_ldm_alloc_stats_t ldm_alloc_stats[SWBWA_CPE_COUNT];
#endif
} swbwa_cpe_task_t;

#endif /* SWBWA_CPE_H */
