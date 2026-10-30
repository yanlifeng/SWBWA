#include "../src/gzip/parallel_pread.h"
#include <cassert>
#include <fcntl.h>
#include <cstdio>

int main() {
    char path[] = "/tmp/swbwa-parallel-pread-XXXXXX";
    const int fd = mkstemp(path);
    assert(fd >= 0);
    std::vector<unsigned char> expected((5U << 20) + 37), actual(expected.size() + 64, 0x91);
    for (size_t i = 0; i < expected.size(); ++i) expected[i] = (i * 37 + i / 8191) & 255;
    for (size_t n = 0; n < expected.size();) {
        const ssize_t wrote = write(fd, expected.data() + n, expected.size() - n);
        assert(wrote > 0); n += size_t(wrote);
    }
    for (unsigned workers : {1U, 2U, 3U, 6U}) {
        SwbwaParallelPread pool(workers);
        for (unsigned iteration = 0; iteration < 32; ++iteration) {
            const size_t offset = iteration * 13, size = expected.size() - offset;
            auto result = pool.read(fd, offset, actual.data() + 32, size);
            assert(result.calls >= workers);
            assert(!std::memcmp(actual.data() + 32, expected.data() + offset, size));
            for (size_t i = 0; i < 32; ++i) assert(actual[i] == 0x91);
            for (size_t i = expected.size() + 32; i < actual.size(); ++i) assert(actual[i] == 0x91);
        }
        bool rejected = false;
        try { pool.read(-1, 0, actual.data(), expected.size()); }
        catch (const std::runtime_error &) { rejected = true; }
        assert(rejected);
        rejected = false;
        try { pool.read(fd, 1, actual.data(), expected.size()); }
        catch (const std::runtime_error &) { rejected = true; }
        assert(rejected);
        pool.read(fd, 0, actual.data(), expected.size());
        assert(actual[0] == expected[0]);
        std::fill(actual.begin(), actual.end(), 0x91);
    }
    close(fd); unlink(path);
    std::puts("PARALLEL_PREAD_PASS workers=1,2,3,6 reuse=32 errors=EBADF,EOF");
}
