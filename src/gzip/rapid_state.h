#ifndef SWBWA_RAPID_STATE_H
#define SWBWA_RAPID_STATE_H
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#if defined(SWBWA_CODEC_MPE) && (!defined(SWBWA_CODEC_LDM) || SWBWA_CODEC_LDM != 0)
#error "SWBWA_CODEC_MPE requires SWBWA_CODEC_LDM=0 to embed distance and history buffers"
#endif
#if defined(SWBWA_CODEC_MPE) && !defined(SWBWA_CODEC_CPE)
#include <cstdlib>
#include <limits>
#endif
#if defined(SWBWA_CODEC_CPE) && SWBWA_CODEC_LDM == 2
extern "C" void *ldm_malloc(size_t);
extern "C" void ldm_free(void *, size_t);
#endif

// Fall back to the existing heap if the hot state would crowd the CPE stack.
template<class T>
class SwbwaRapidState {
public:
    SwbwaRapidState() {
#if defined(SWBWA_CODEC_MPE) && !defined(SWBWA_CODEC_CPE)
        static_assert(sizeof(T) <= std::numeric_limits<size_t>::max() - (alignof(T) - 1),
                      "Aligned state allocation size overflows");
        size_t space = sizeof(T) + alignof(T) - 1;
        // Native aligned new/delete can retain MPE memory with xmalloc enabled.
        // Pair malloc/free explicitly and keep ownership during construction.
        std::unique_ptr<void, decltype(&std::free)> allocation(std::malloc(space), &std::free);
        if (!allocation) throw std::bad_alloc();
        void *aligned = allocation.get();
        if (!std::align(alignof(T), sizeof(T), aligned, space)) throw std::bad_alloc();
        value_ = ::new (aligned) T;
        backing_ = allocation.release();
#else
#if defined(SWBWA_CODEC_CPE) && SWBWA_CODEC_LDM == 2
        backing_ = ldm_malloc(sizeof(T) + 63);
        volatile unsigned char marker;
        uintptr_t start = reinterpret_cast<uintptr_t>(backing_);
        uintptr_t stack = reinterpret_cast<uintptr_t>(&marker);
        if (backing_ && (start >= stack || sizeof(T) + 63 > stack - start ||
            stack - start - sizeof(T) - 63 < (72U << 10))) {
            ldm_free(backing_, sizeof(T) + 63); backing_ = nullptr;
        }
        if (backing_) {
            try { value_ = new (reinterpret_cast<void *>((start + 63) & ~uintptr_t(63))) T; }
            catch (...) { ldm_free(backing_, sizeof(T) + 63); throw; }
        } else
#endif
            heap_.reset(new T);
        if (!value_) value_ = heap_.get();
#endif
    }
    ~SwbwaRapidState() {
#if defined(SWBWA_CODEC_MPE) && !defined(SWBWA_CODEC_CPE)
        value_->~T();
        std::free(backing_);
#else
#if defined(SWBWA_CODEC_CPE) && SWBWA_CODEC_LDM == 2
        if (backing_) { value_->~T(); ldm_free(backing_, sizeof(T) + 63); }
#endif
#endif
    }
    SwbwaRapidState(const SwbwaRapidState &) = delete;
    SwbwaRapidState &operator=(const SwbwaRapidState &) = delete;
    T *operator->() { return value_; }
    bool in_ldm() const {
#ifdef SWBWA_CODEC_CPE
        return reinterpret_cast<uintptr_t>(value_) < (256U << 10);
#else
        return false;
#endif
    }
private:
#if !defined(SWBWA_CODEC_MPE) || defined(SWBWA_CODEC_CPE)
    std::unique_ptr<T> heap_;
#endif
    void *backing_ = nullptr;
    T *value_ = nullptr;
};
#endif
