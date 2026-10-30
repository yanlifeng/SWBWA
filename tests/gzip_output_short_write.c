#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/uio.h>
#include <unistd.h>

ssize_t __real_writev(int, const struct iovec *, int);
ssize_t __real_pwritev64(int, const struct iovec *, int, off_t);
long __real_sysconf(int);

long __wrap_sysconf(int name)
{
    const char *option = getenv("SWBWA_TEST_IOV_LIMIT");
    return name == _SC_IOV_MAX && option ? strtol(option, NULL, 10) : __real_sysconf(name);
}

static int limited(const struct iovec *input, int count, struct iovec *output)
{
    const char *option = getenv("SWBWA_TEST_WRITE_LIMIT");
    size_t remaining = option ? (size_t)strtoul(option, NULL, 10) : (size_t)-1;
    int i, used = 0;
    static __thread int interrupted;
    long limit = __wrap_sysconf(_SC_IOV_MAX);
    assert(count > 0 && count <= 1024 && count <= (limit > 0 ? limit : 16));
    if (!interrupted && getenv("SWBWA_TEST_WRITE_EINTR")) {
        interrupted = 1;
        errno = EINTR;
        return -1;
    }
    for (i = 0; i < count && remaining; ++i) {
        output[used] = input[i];
        if (output[used].iov_len > remaining) output[used].iov_len = remaining;
        remaining -= output[used].iov_len;
        ++used;
    }
    return used;
}

ssize_t __wrap_writev(int fd, const struct iovec *vectors, int count)
{
    struct iovec parts[1024];
    int n = limited(vectors, count, parts);
    return n < 0 ? -1 : __real_writev(fd, parts, n);
}

ssize_t __wrap_pwritev64(int fd, const struct iovec *vectors, int count, off_t offset)
{
    struct iovec parts[1024];
    int n = limited(vectors, count, parts);
    return n < 0 ? -1 : __real_pwritev64(fd, parts, n, offset);
}
