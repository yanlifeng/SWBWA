#ifndef SWBWA_MPI_H
#define SWBWA_MPI_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

#include "swbwa_config.h"

typedef struct {
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
void swbwa_mpi_print_rank_ordered(void (*printer)(void));
void swbwa_mpi_abort(const char *message);

int swbwa_mpi_fastq_range(const char *read1_path, const char *read2_path,
                          swbwa_fastq_range_t *range);

int swbwa_fastq_estimate_chunk_bytes(const char *read1_path,
                                     const char *read2_path,
                                     int64_t target_bases,
                                     int64_t *file_size,
                                     int64_t *chunk_bytes);

#if SWBWA_ENABLE_MPI_FASTQ_SCHEDULER
int swbwa_mpi_fastq_scheduler_open(const char *read1_path,
                                   const char *read2_path,
                                   int64_t target_bases,
                                   swbwa_fastq_range_t *assigned_range,
                                   int64_t *chunk_count);
int64_t swbwa_mpi_fastq_scheduler_chunk_bytes(void);
int swbwa_mpi_fastq_scheduler_next(swbwa_fastq_range_t *range);
void swbwa_mpi_fastq_scheduler_add_records(int64_t records);
void swbwa_mpi_fastq_scheduler_stats(int64_t *chunks, int64_t *records,
                                     int64_t *bytes);
void swbwa_mpi_fastq_scheduler_close(void);
#endif

int swbwa_input_register_fd(int fd, int64_t start, int64_t end);
int swbwa_input_set_range(int fd, int64_t start, int64_t end);
void swbwa_input_unregister_fd(int fd);
ssize_t swbwa_input_read(int fd, void *buffer, size_t bytes);

void swbwa_mpi_call_lock(void);
void swbwa_mpi_call_unlock(void);

#endif /* SWBWA_MPI_H */
