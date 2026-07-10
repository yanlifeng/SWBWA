#ifndef SWBWA_CPE_H
#define SWBWA_CPE_H

#include "swbwa_config.h"

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

    int *completion_flags;
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
} swbwa_cpe_task_t;

#endif /* SWBWA_CPE_H */
