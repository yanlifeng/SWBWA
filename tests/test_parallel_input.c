#include "swbwa_input.h"
#include "swbwa_host_workers.h"
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef TEST_WRAP_READ
static int read_mode;
static _Thread_local unsigned int read_attempts;
ssize_t __real_pread(int fd, void *buffer, size_t count, off_t offset);
ssize_t __wrap_pread(int fd, void *buffer, size_t count, off_t offset)
{
    if (read_mode == 2) { errno = EIO; return -1; }
    if (read_mode == 1) {
        if (++read_attempts % 11 == 1) { errno = EINTR; return -1; }
        if (count > 7) count = 7;
    }
    return __real_pread(fd, buffer, count, offset);
}
#endif

static void check_read(int fd, int readers, const unsigned char *source,
                       size_t first, size_t second, int profile)
{
    unsigned char buffers[2][8192];
    swbwa_input_slice_t inputs[2] = {
        {fd, (char *)buffers[0] + 16, 3, first},
        {fd, (char *)buffers[1] + 16, 129, second}
    };
    swbwa_input_result_t result;
    int i;
    memset(buffers, 0xa5, sizeof(buffers));
    assert(lseek(fd, 29, SEEK_SET) == 29);
    assert(swbwa_input_read(inputs, 2, readers, profile, &result) == 0);
    assert(lseek(fd, 0, SEEK_CUR) == 29);
    assert(result.bytes == first + second);
    for (i = 0; i < 2; ++i) {
        size_t j;
        assert(memcmp(inputs[i].buffer, source + inputs[i].offset, inputs[i].length) == 0);
        for (j = 0; j < 16; ++j) assert(buffers[i][j] == 0xa5);
        for (j = 16 + inputs[i].length; j < sizeof(buffers[i]); ++j)
            assert(buffers[i][j] == 0xa5);
    }
}

int main(int argc, char **argv)
{
    const char *invalid_readers[] = {"", "-1", "+1", "7", "16", " 6", "6 "};
    unsigned char source[16384];
    swbwa_input_slice_t input;
    swbwa_input_result_t result;
    int fd, readers, i;
    char byte;
    assert(argc == 2);
    assert(swbwa_input_reader_count(NULL) == (SWBWA_HOST_MPE_THREADS > 1 ? 6 : 0));
    for (i = 0; i <= 6; ++i) {
        char value[2] = {(char)('0' + i), '\0'};
        if (i <= SWBWA_HOST_MPE_THREADS)
            assert(swbwa_input_reader_count(value) == i);
        else
            assert(swbwa_input_reader_count(value) == -1 && errno == EINVAL);
    }
    for (i = 0; i < (int)(sizeof(invalid_readers) / sizeof(invalid_readers[0])); ++i)
        assert(swbwa_input_reader_count(invalid_readers[i]) == -1 && errno == EINVAL);
    fd = open(argv[1], O_RDWR | O_CREAT | O_EXCL, 0600);
    assert(fd >= 0);
    for (i = 0; i < (int)sizeof(source); ++i) source[i] = (unsigned char)(i * 31 + i / 17);
    assert(write(fd, source, sizeof(source)) == sizeof(source));
    assert(swbwa_host_workers_init() == 0);
    readers = swbwa_input_reader_count(NULL);
    if (readers > 0) check_read(fd, readers, source, 4097, 3001, 1);
    for (readers = 1; readers <= SWBWA_HOST_MPE_THREADS; ++readers) {
        check_read(fd, readers, source, 1, 2, 1);
        check_read(fd, readers, source, 4097, 3001, 1);
        check_read(fd, readers, source, 4096, 0, 0);
        check_read(fd, readers, source, 0, 0, 1);
#ifdef TEST_WRAP_READ
        read_mode = 1;
        check_read(fd, readers, source, 1001, 557, 1);
        read_mode = 0;
#endif
    }
    input.fd = fd;
    input.buffer = &byte;
    input.offset = (INT64_C(1) << 32) + 17;
    input.length = 1;
    assert(pwrite(fd, "x", 1, (off_t)input.offset) == 1);
    assert(swbwa_input_read(&input, 1, 1, 1, &result) == 0 && byte == 'x');
    ++input.offset;
    assert(swbwa_input_read(&input, 1, 1, 1, &result) == -1 && errno == EIO);
    input.offset = INT64_MAX;
    assert(swbwa_input_read(&input, 1, 1, 1, &result) == -1 && errno == EOVERFLOW);
    input.offset = 0;
    assert(swbwa_input_read(&input, 1, 0, 1, &result) == -1 && errno == EINVAL);
    assert(swbwa_input_read(&input, 1, SWBWA_HOST_MPE_THREADS + 1, 1, &result) == -1);
#ifdef TEST_WRAP_READ
    read_mode = 2;
    assert(swbwa_input_read(&input, 1, SWBWA_HOST_MPE_THREADS, 1, &result) == -1 && errno == EIO);
    read_mode = 0;
#endif
    input.length = 0;
    input.buffer = NULL;
    input.fd = -1;
    assert(swbwa_input_read(&input, 1, 1, 1, &result) == 0);
    swbwa_host_workers_destroy();
    assert(close(fd) == 0);
    assert(unlink(argv[1]) == 0);
    puts("parallel positioned input: PASS");
    return 0;
}
