#ifndef SWBWA_RAPID_READER_H
#define SWBWA_RAPID_READER_H
#include <cstddef>
#include <memory>
#include <vector>
#include "decoded_view.h"
struct SwbwaGzipMetrics;
class SwbwaLineIndex;

class SwbwaRapidReader {
public:
    SwbwaRapidReader(const char *path, size_t chunk, unsigned workers, SwbwaGzipMetrics &metrics,
                    const std::vector<unsigned char> *memory = nullptr, bool views = false);
    ~SwbwaRapidReader();
    size_t read(void *out, size_t bytes, SwbwaLineIndex *lines = nullptr);
    bool eof() const;
    bool peek(size_t maximum, SwbwaDecodedView &view);
    void consume(size_t bytes);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
#endif
