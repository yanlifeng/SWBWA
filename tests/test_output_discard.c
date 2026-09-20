#include "swbwa_output.h"

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int swbwa_mpi_rank(void) { return 0; }
int swbwa_mpi_size(void) { return 1; }
int swbwa_mpi_is_root(void) { return 1; }
void swbwa_mpi_print_rank_ordered(void (*printer)(void)) { printer(); }

int main(int argc, char **argv)
{
    static const char sam[] = "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    unsigned char unaligned[1026];
    size_t i;

    assert(argc == 2);
    assert(access(argv[1], F_OK) != 0);
    assert(swbwa_output_write(sam, sizeof(sam) - 1) == -1);
    assert(swbwa_output_open(argv[1], 1) == 0);
    assert(strcmp(swbwa_output_name(), "(discard)") == 0);
    assert(swbwa_output_open(argv[1], 1) == -1 && errno == EALREADY);
    assert(swbwa_output_write(NULL, 0) == 0);
    assert(swbwa_output_write(NULL, 1) == -1 && errno == EINVAL);
    assert(swbwa_output_write(sam, sizeof(sam) - 1) == 0);
    for (i = 0; i < sizeof(unaligned); ++i) unaligned[i] = (unsigned char)i;
    for (i = 1; i <= 1024; ++i)
        assert(swbwa_output_write(unaligned + 1, i) == 0);
    assert(swbwa_output_flush() == 0);
    assert(swbwa_output_close() == 0);
    assert(swbwa_output_close() == 0);
    assert(access(argv[1], F_OK) != 0);
    return 0;
}
