#ifndef SWBWA_MPI_DECODER_PROTOCOL_H
#define SWBWA_MPI_DECODER_PROTOCOL_H

// Derived from SWBWA's distributed gzip protocol experiment. Generalized to
// rank-local job queues; the underlying codec retains its own provenance.
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace swbwa_mpi_gzip {
using std::size_t;
constexpr size_t HISTORY = 32768;
constexpr size_t PAGE = 4096;
enum Field { PRESENT, OK, START, START_MAX, END, BYTES, FINAL, OUTPUT_FULL,
             BASE, LIMIT, KNOWN, FIELDS };
using Meta = std::array<uint64_t, FIELDS>;
using Window = std::vector<unsigned char>;
using Plan = std::array<uint64_t, 3>; // selected, global offset, true start bit
static_assert(sizeof(Meta) == FIELDS * sizeof(uint64_t), "wire metadata must be contiguous");
static_assert(sizeof(Plan) == 3 * sizeof(uint64_t), "wire plan must be contiguous");

inline uint64_t add(uint64_t a, uint64_t b) {
    if (b > std::numeric_limits<uint64_t>::max() - a)
        throw std::runtime_error("integer overflow");
    return a + b;
}
inline uint64_t multiply(uint64_t a, uint64_t b) {
    if (b && a > std::numeric_limits<uint64_t>::max() / b)
        throw std::runtime_error("integer overflow");
    return a * b;
}

struct Chain {
    std::vector<size_t> jobs;
    uint64_t cursor, decoded;
    bool final;
};

inline Chain select(const std::vector<Meta> &metas, std::vector<Plan> &plan,
                    unsigned ranks, unsigned count, unsigned anchor,
                    uint64_t cursor, uint64_t decoded, uint64_t file_bits,
                    uint64_t capacity) {
    if (!ranks || !count || anchor >= ranks || metas.size() != size_t(ranks) * count)
        throw std::runtime_error("invalid job grid");
    plan.assign(metas.size(), Plan{});
    for (const auto &m : metas) {
        if (!m[PRESENT] || !m[OK]) continue;
        if (m[BASE] > m[START] || m[START] > m[START_MAX] ||
            m[START_MAX] >= m[END] || m[END] > m[LIMIT] ||
            m[LIMIT] > file_bits || m[BYTES] > capacity || m[FINAL] > 1)
            throw std::runtime_error("invalid candidate metadata");
    }
    const Meta &first = metas[size_t(anchor) * count];
    if (!first[PRESENT] || !first[OK])
        throw std::runtime_error(first[OUTPUT_FULL] ? "anchored DISCOVER exhausted memory budget"
            : "anchored DISCOVER failed: corrupt stream or insufficient lookahead");
    if (!first[KNOWN] || first[START] != cursor || first[START_MAX] != cursor)
        throw std::runtime_error("anchor does not match confirmed cursor");
    Chain chain{{}, cursor, decoded, false};
    // Jobs occupy monotonically increasing, disjoint search intervals. A long
    // anchored block may skip many intervals; only exact continuity commits.
    for (unsigned group = 0; group < ranks; ++group) {
        const unsigned rank = (anchor + group) % ranks;
        for (unsigned j = 0; j < count; ++j) {
            const size_t id = size_t(rank) * count + j;
            const auto &m = metas[id];
            if (!m[PRESENT] || !m[OK] || chain.cursor < m[START] || chain.cursor > m[START_MAX])
                continue;
            plan[id] = {1, chain.decoded, chain.cursor};
            chain.jobs.push_back(id);
            chain.decoded = add(chain.decoded, m[BYTES]);
            chain.cursor = m[END];
            chain.final = m[FINAL] != 0;
            if (chain.final) return chain;
        }
    }
    if (chain.jobs.empty() || chain.cursor <= cursor)
        throw std::runtime_error("candidate chain made no bit progress");
    return chain;
}

inline Window advance(const Window &previous, const uint16_t *symbols, size_t size) {
    if (previous.size() > HISTORY) throw std::runtime_error("oversized dictionary");
    const size_t n = std::min(size, HISTORY);
    const size_t keep = size < HISTORY ? std::min(previous.size(), HISTORY - n) : 0;
    Window next;
    next.reserve(keep + n);
    next.insert(next.end(), previous.end() - keep, previous.end());
    for (size_t i = size - n; i < size; ++i) {
        unsigned symbol = symbols[i];
        if (symbol >= 256) {
            if (symbol < 65536U - previous.size())
                throw std::runtime_error("invalid tail marker / missing history");
            symbol = previous[symbol - (65536U - previous.size())];
        }
        next.push_back(static_cast<unsigned char>(symbol));
    }
    return next;
}
} // namespace swbwa_mpi_gzip
#endif
