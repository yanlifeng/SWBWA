#ifndef SWBWA_MPI_H
#define SWBWA_MPI_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

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
void swbwa_mpi_print_rank_ordered(void (*printer)(void));
void swbwa_mpi_abort(const char *message);

int swbwa_mpi_fastq_range(const char *read1_path, const char *read2_path,
                          swbwa_fastq_range_t *range);

int swbwa_input_register_fd(int fd, int64_t start, int64_t end);
void swbwa_input_unregister_fd(int fd);
ssize_t swbwa_input_read(int fd, void *buffer, size_t bytes);

#endif /* SWBWA_MPI_H */
