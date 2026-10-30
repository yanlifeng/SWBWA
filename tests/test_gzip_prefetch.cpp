// Exercise the private Source without adding a production testing API.
#include "../src/gzip/mpi_decoder.cpp"
#include <cassert>
#include <dirent.h>
#include <set>
#include <sys/syscall.h>

static const std::thread::id owner = std::this_thread::get_id();
static std::mutex observation_mutex;
static std::set<long> worker_ids;
static std::atomic<unsigned> active_reads{0}, completed_reads{0};
static std::atomic<bool> slow_read{false};

extern "C" ssize_t __real_pread64(int, void *, size_t, off64_t);
extern "C" ssize_t __wrap_pread64(int fd, void *buffer, size_t size, off64_t offset) {
    const bool worker = std::this_thread::get_id() != owner;
    if (worker) {
        std::lock_guard<std::mutex> lock(observation_mutex);
        worker_ids.insert(syscall(SYS_gettid));
        ++active_reads;
    }
    if (worker && slow_read) usleep(10000);
    const ssize_t result = __real_pread64(fd, buffer, size, offset);
    const int error = errno;
    if (worker) {
        assert(result >= 0); // In particular, close must not race the read.
        ++completed_reads;
        --active_reads;
    }
    errno = error;
    return result;
}

extern "C" double __wrap_MPI_Wtime() {
    assert(std::this_thread::get_id() == owner);
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

extern "C" void swbwa_codec_launch(swbwa_codec_batch *) { std::abort(); }

static size_t worker_count() {
    std::lock_guard<std::mutex> lock(observation_mutex);
    return worker_ids.size();
}

static size_t thread_count() {
    DIR *directory = opendir("/proc/self/task");
    assert(directory);
    size_t count = 0;
    while (dirent *entry = readdir(directory)) if (entry->d_name[0] != '.') ++count;
    closedir(directory);
    return count;
}

static void joined(size_t baseline) {
    // The kernel can briefly retain a task entry after pthread_join returns.
    for (unsigned i = 0; i < 1000 && thread_count() != baseline; ++i) usleep(1000);
    assert(thread_count() == baseline && active_reads == 0);
}

template<class Exception, class Function> static void rejects(Function function) {
    bool rejected = false;
    try { function(); } catch (const Exception &) { rejected = true; }
    assert(rejected);
}

static void check(Source &source, const std::vector<unsigned char> &expected,
                  size_t offset, size_t size) {
    auto bytes = source.read_batch(offset, size);
    assert(!std::memcmp(bytes.get(), expected.data() + offset, size));
    for (size_t i = 0; i < GUARD; ++i) assert(bytes[size + i] == 0);
}

int main() {
    char path[] = "/tmp/swbwa-prefetch-XXXXXX";
    const int fd = mkstemp(path);
    assert(fd >= 0);
    std::vector<unsigned char> expected(2U << 20);
    for (size_t i = 0; i < expected.size(); ++i) expected[i] = (i * 17 + i / 31) & 255;
    size_t written = 0;
    while (written < expected.size()) {
        const ssize_t n = write(fd, expected.data() + written, expected.size() - written);
        assert(n > 0);
        written += size_t(n);
    }
    close(fd);
    const size_t baseline = thread_count();
    assert(setenv("SWBWA_GZIP_PREFETCH_MIB", "1", 1) == 0);
    {
        swbwa_mpi_decoder_metrics a{}, b{};
        Source first(path, a), second(path, b);
        first.prefetch(0, 0);
        assert(thread_count() == baseline && worker_count() == 0);
        for (size_t round = 0; round < 512; ++round) {
            const size_t offset = round * 37, size = 4096 + round;
            first.prefetch(offset, size);
            second.prefetch(offset + 1, size);
            check(first, expected, offset, size);
            check(second, expected, offset + 1, size);
            assert(worker_count() == 2);
        }
        assert(thread_count() == baseline + 2);
        assert(a.prefetch_batches == 512 && b.prefetch_batches == 512);
        first.prefetch(17, 4096);
        rejects<std::logic_error>([&] { first.prefetch(17, 4096); });
        rejects<std::logic_error>([&] { first.read_batch(18, 4096); });
        rejects<std::logic_error>([&] { first.read_batch(17, 4095); });
        check(first, expected, 17, 4096);
        const uint64_t before = a.prefetch_bytes;
        first.prefetch(0, expected.size());
        check(first, expected, 0, expected.size());
        assert(a.prefetch_bytes - before == (1U << 20));
        first.prefetch(expected.size() + 1, 64);
        rejects<std::runtime_error>([&] { first.read_batch(expected.size() + 1, 64); });
        first.prefetch(0, 64);
        check(first, expected, 0, 64);
        assert(worker_count() == 2);
    }
    joined(baseline);
    {
        swbwa_mpi_decoder_metrics metrics{};
        Source source(path, metrics);
        struct Range { size_t old_offset, old_size, offset, size, cached; };
        const Range ranges[] = {
            {0, 4096, 512, 4096, 3584},
            {512, 4096, 0, 4096, 3584},
            {0, 4096, 0, 4096, 4096},
            {0, 4096, 8192, 4096, 0},
            {1024, 1024, 0, 4096, 1024},
            {0, expected.size(), 0, expected.size(), 1U << 20},
        };
        for (const auto &range : ranges) {
            std::vector<unsigned char> old(expected.begin() + range.old_offset,
                                            expected.begin() + range.old_offset + range.old_size);
            const uint64_t read_before = metrics.source_read_bytes;
            const uint64_t cached_before = metrics.source_cache_bytes;
            source.prefetch(range.offset, range.size, old.data(), range.old_offset, old.size());
            // The asynchronous worker must not retain a pointer to the old window.
            std::fill(old.begin(), old.end(), 0);
            std::vector<unsigned char>().swap(old);
            check(source, expected, range.offset, range.size);
            assert(metrics.source_read_bytes - read_before == range.size - range.cached);
            assert(metrics.source_cache_bytes - cached_before == range.cached);
        }
        assert(setenv("SWBWA_GZIP_REUSE_INPUT", "0", 1) == 0);
        const uint64_t read_before = metrics.source_read_bytes;
        const uint64_t cached_before = metrics.source_cache_bytes;
        source.prefetch(0, 4096, expected.data(), 0, expected.size());
        check(source, expected, 0, 4096);
        assert(metrics.source_read_bytes - read_before == 4096);
        assert(metrics.source_cache_bytes == cached_before);
        assert(unsetenv("SWBWA_GZIP_REUSE_INPUT") == 0);
    }
    joined(baseline);
    assert(setenv("SWBWA_GZIP_PREFETCH_MIB", "0", 1) == 0);
    {
        swbwa_mpi_decoder_metrics metrics{};
        Source source(path, metrics);
        source.prefetch(0, 64);
        check(source, expected, 0, 64);
        assert(metrics.prefetch_batches == 0);
        assert(thread_count() == baseline);
    }
    assert(setenv("SWBWA_GZIP_PREFETCH_MIB", "1", 1) == 0);
    slow_read = true;
    for (unsigned round = 0; round < 32; ++round) {
        const unsigned before = completed_reads;
        {
            swbwa_mpi_decoder_metrics metrics{};
            Source source(path, metrics);
            // Close after an abandoned result, abandoned error, or consumed error.
            const size_t offset = round % 3 ? expected.size() + 1 : 0;
            source.prefetch(offset, 64);
            if (round % 3 == 2)
                rejects<std::runtime_error>([&] { source.read_batch(offset, 64); });
        }
        assert(completed_reads == before + 1);
        joined(baseline);
    }
    assert(unlink(path) == 0);
    std::puts("GZIP_PREFETCH_PASS: 2 workers / 1024 requests; 32 create-close cycles; range, prefix, reuse, error, disabled");
}
