// Test-only compressed-byte replay cache. A full warm-up decode captures only
// this rank's pread requests; timed passes reject misses instead of touching disk.
// Link with --wrap=pread and --wrap=pread64. Not part of the application build.
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

extern "C" ssize_t __real_pread(int, void *, size_t, off_t);
extern "C" ssize_t __real_pread64(int, void *, size_t, off64_t);
namespace {
std::mutex mutex;
using Key = std::pair<uint64_t, size_t>;
std::map<Key, std::vector<unsigned char>> blocks;
struct stat source;
unsigned mode;
uint64_t bytes, hits, misses;
uint64_t limit = UINT64_C(2) << 30;
}

extern "C" int swbwa_benchmark_cache_open(const char *path) {
    if (stat(path, &source)) return -1;
    const char *setting = std::getenv("SWBWA_HISTORY_CACHE_MIB");
    if (setting) {
        char *end = nullptr;
        const unsigned long mib = std::strtoul(setting, &end, 10);
        if (!*setting || *end || mib < 1 || mib > 8192) { errno = EINVAL; return -1; }
        limit = uint64_t(mib) << 20;
    }
    mode = 1;
    return 0;
}

extern "C" void swbwa_benchmark_cache_freeze(int rank) {
    std::lock_guard<std::mutex> guard(mutex);
    std::printf("HISTORY_CACHE rank=%d bytes=%llu blocks=%zu warmup_hits=%llu warmup_misses=%llu\n",
        rank, (unsigned long long)bytes, blocks.size(), (unsigned long long)hits, (unsigned long long)misses);
    mode = 2;
}

static ssize_t cached_read(int fd, void *buffer, size_t count, off64_t offset) {
    if (!mode) return __real_pread64(fd, buffer, count, offset);
    std::unique_lock<std::mutex> guard(mutex);
    struct stat info;
    if (fstat(fd, &info)) return -1;
    if (info.st_dev != source.st_dev || info.st_ino != source.st_ino)
        return __real_pread64(fd, buffer, count, offset);
    const Key key{uint64_t(offset), count};
    auto found = blocks.find(key);
    if (found == blocks.end()) {
        ++misses;
        if (mode == 2 || count > limit - bytes) {
            std::fprintf(stderr, "HISTORY_CACHE_ERROR miss/limit offset=%lld count=%zu mode=%u\n",
                (long long)offset, count, mode);
            errno = EIO; return -1;
        }
        std::vector<unsigned char> data(count);
        const ssize_t n = __real_pread64(fd, data.data(), count, offset);
        if (n < 0) return n;
        data.resize(size_t(n)); bytes += data.size();
        found = blocks.emplace(key, std::move(data)).first;
    } else ++hits;
    const auto &data = found->second;
    // Entries never change or move after insertion. Do not serialize the
    // application's parallel copies with this test-only lookup mutex.
    guard.unlock();
    if (!data.empty()) __builtin_memcpy(buffer, data.data(), data.size());
    return ssize_t(data.size());
}

extern "C" ssize_t __wrap_pread(int fd, void *buffer, size_t count, off_t offset) {
    return cached_read(fd, buffer, count, offset);
}
extern "C" ssize_t __wrap_pread64(int fd, void *buffer, size_t count, off64_t offset) {
    return cached_read(fd, buffer, count, offset);
}
