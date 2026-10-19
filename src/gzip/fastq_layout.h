#ifndef SWBWA_GZIP_FASTQ_LAYOUT_H
#define SWBWA_GZIP_FASTQ_LAYOUT_H

#include "line_index.h"
#include <string>

/* Match the legacy CPE formatter's byte partition, including its lookbehind.
 * Equal mate counts alone do not make shared byte offsets safe. */
inline size_t swbwa_fastq_legacy_boundary(const unsigned char *data,
                                         const SwbwaLineIndex &index,
                                         size_t position)
{
    const size_t size = index.bytes();
    size_t line = index.lines_before(data, position);
    auto skip = [&]() {
        position = line < index.lines() ? index.end_of_line(data, ++line) : size;
    };
    skip();
    unsigned searched = 0;
    while (position < size && data[position] != '@') {
        if (++searched > 4)
            throw std::runtime_error("invalid FASTQ near CPE partition boundary");
        skip();
    }
    if (position == size) return size;
    const size_t candidate = position;
    skip();
    if (position == size) return candidate;
    if (data[position] == '@') return position;
    skip();
    if (position >= size || data[position] != '+')
        throw std::runtime_error("invalid FASTQ near CPE partition boundary");
    return candidate;
}

inline void swbwa_fastq_check_layout(const unsigned char *r1,
                                     const SwbwaLineIndex &a,
                                     const unsigned char *r2,
                                     const SwbwaLineIndex &b, unsigned count)
{
    if (!count || count > 384)
        throw std::invalid_argument("CPE count must be 1..384");
    if (!r1 || !r2 || a.bytes() != b.bytes())
        throw std::runtime_error("PE layout incompatible with shared CPE offsets: unequal batch bytes");
    if (a.lines() != b.lines() || a.lines() % 4)
        throw std::runtime_error("PE layout has unequal or incomplete record counts");
    if (a.bytes() < size_t(count) * 1024) return;
    const size_t step = (a.bytes() + count - 1) / count;
    for (unsigned i = 1; i < count; ++i) {
        const size_t nominal = std::min(a.bytes(), size_t(i) * step);
        const size_t boundary = swbwa_fastq_legacy_boundary(r1, a, nominal - 1024);
        const size_t rank = a.lines_before(r1, boundary);
        if (rank % 4 || (boundary && r1[boundary - 1] != '\n') ||
            (rank ? b.end_of_line(r2, rank) : 0) != boundary)
            throw std::runtime_error("PE layout incompatible with shared CPE offsets at partition " +
                                     std::to_string(i));
    }
}

#endif
