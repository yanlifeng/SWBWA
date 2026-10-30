#ifndef SWBWA_PARALLEL_PREAD_H
#define SWBWA_PARALLEL_PREAD_H

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
#include <athread.h>
#endif

// Dedicated I/O workers must not hold the compute-helper dispatch lock while
// blocked in pread. One caller owns each pool; workers never enter MPI/CPE.
class SwbwaParallelPread {
public:
    struct Result { uint64_t calls = 0; double seconds = 0; };
    explicit SwbwaParallelPread(unsigned count) : count_(count), results_(count) {
        if (!count || count > 6) throw std::invalid_argument("invalid pread worker count");
    }
    ~SwbwaParallelPread() { stop(); }
    SwbwaParallelPread(const SwbwaParallelPread &) = delete;
    SwbwaParallelPread &operator=(const SwbwaParallelPread &) = delete;

    Result read(int fd, uint64_t offset, void *buffer, size_t size) {
        if (count_ == 1 || size < (1U << 20)) return serial(fd, offset, buffer, size);
        start();
        const auto begin = std::chrono::steady_clock::now();
        std::unique_lock<std::mutex> lock(mutex_);
        fd_ = fd; offset_ = offset; buffer_ = static_cast<unsigned char *>(buffer); size_ = size;
        error_ = nullptr; remaining_ = count_; ++generation_;
        work_.notify_all();
        done_.wait(lock, [this] { return remaining_ == 0; });
        if (error_) std::rethrow_exception(error_);
        Result result;
        for (const auto &part : results_) result.calls += part.calls;
        result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        return result;
    }

    static Result serial(int fd, uint64_t offset, void *buffer, size_t size) {
        const auto begin = std::chrono::steady_clock::now();
        Result result;
        auto *p = static_cast<unsigned char *>(buffer);
        while (size) {
            ++result.calls;
            const ssize_t n = pread(fd, p, size, static_cast<off_t>(offset));
            if (n < 0 && errno == EINTR) continue;
            if (n < 0) throw std::runtime_error(std::string("pread gzip: ") + std::strerror(errno));
            if (!n) throw std::runtime_error("truncated or changed input");
            size -= size_t(n); p += n; offset += uint64_t(n);
        }
        result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        return result;
    }

private:
    void start() {
        if (!threads_.empty()) return;
        try {
            for (unsigned i = 0; i < count_; ++i) threads_.emplace_back(&SwbwaParallelPread::worker, this, i);
        } catch (...) { stop(); throw; }
    }
    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
            work_.notify_all();
        }
        for (auto &thread : threads_) if (thread.joinable()) thread.join();
    }
    void worker(unsigned id) {
        std::exception_ptr binding_error;
#if defined(__sw_host__) && SWBWA_HOST_MPE_THREADS > 1
        try {
            SET_TO_SWCG(id);
            unsigned long actual;
            __asm__ volatile("rcid %0" : "=r"(actual));
            if ((actual & 7) != id) throw std::runtime_error("gzip pread worker MPE binding failed");
        } catch (...) { binding_error = std::current_exception(); }
#endif
        unsigned long generation = 0;
        std::unique_lock<std::mutex> lock(mutex_);
        for (;;) {
            work_.wait(lock, [this, generation] { return stopping_ || generation_ != generation; });
            if (stopping_) return;
            generation = generation_;
            const size_t begin = size_ / count_ * id + std::min<size_t>(size_ % count_, id);
            const size_t end = size_ / count_ * (id + 1) + std::min<size_t>(size_ % count_, id + 1);
            lock.unlock();
            Result result;
            std::exception_ptr error = binding_error;
            if (!error) try { result = serial(fd_, offset_ + begin, buffer_ + begin, end - begin); }
            catch (...) { error = std::current_exception(); }
            lock.lock();
            results_[id] = result;
            if (error && !error_) error_ = error;
            if (--remaining_ == 0) done_.notify_one();
        }
    }
    unsigned count_, remaining_ = 0;
    unsigned long generation_ = 0;
    bool stopping_ = false;
    int fd_ = -1;
    uint64_t offset_ = 0;
    unsigned char *buffer_ = nullptr;
    size_t size_ = 0;
    std::mutex mutex_;
    std::condition_variable work_, done_;
    std::exception_ptr error_;
    std::vector<Result> results_;
    std::vector<std::thread> threads_;
};
#endif
