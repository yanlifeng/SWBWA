#ifndef SWBWA_MPI_H
#define SWBWA_MPI_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

#include "swbwa_config.h"

typedef struct {
    int64_t chunk_id;
    int64_t start;
    int64_t end;
    int64_t file_size;
    int64_t first_record;
    int64_t record_count;
} swbwa_fastq_range_t;

int swbwa_mpi_init(int *argc, char ***argv);
void swbwa_mpi_finalize(void);
int swbwa_mpi_rank(void);
int swbwa_mpi_size(void);
int swbwa_mpi_is_root(void);
int swbwa_mpi_barrier(void);
int swbwa_mpi_progress(void);
int swbwa_mpi_progress_thread_start(int debug_enabled);
int swbwa_mpi_progress_thread_active(void);
int swbwa_mpi_progress_thread_stop(void);
void swbwa_mpi_progress_thread_report(void);
void swbwa_mpi_print_rank_ordered(void (*printer)(void));
void swbwa_mpi_abort(const char *message);

int swbwa_mpi_fastq_range(const char *read1_path, const char *read2_path,
                          swbwa_fastq_range_t *range);

int swbwa_fastq_chunk_bytes(const char *read1_path,
                            const char *read2_path,
                            int64_t bytes_per_cg,
                            int64_t *file_size,
                            int64_t *chunk_bytes);

#if SWBWA_USE_MPI && \
    SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC
int swbwa_mpi_fastq_scheduler_open(const char *read1_path,
                                   const char *read2_path,
                                   int64_t bytes_per_cg,
                                   int debug_enabled,
                                   swbwa_fastq_range_t *assigned_range,
                                   int64_t *chunk_count);
int64_t swbwa_mpi_fastq_scheduler_chunk_bytes(void);
int swbwa_mpi_fastq_scheduler_tail_percent(void);
int64_t swbwa_mpi_fastq_scheduler_micro_chunk_bytes(void);
int64_t swbwa_mpi_fastq_scheduler_fine_chunk_bytes(void);
int swbwa_mpi_fastq_scheduler_fine_tail_waves(void);
int swbwa_mpi_fastq_scheduler_next(swbwa_fastq_range_t *range);
void swbwa_mpi_fastq_scheduler_record_stage2(int64_t chunk_id,
                                             int64_t records,
                                             double seconds);
void swbwa_mpi_fastq_scheduler_stats(int64_t *chunks, int64_t *records,
                                     int64_t *bytes);
void swbwa_mpi_fastq_scheduler_close(void);
#endif

#endif /* SWBWA_MPI_H */
