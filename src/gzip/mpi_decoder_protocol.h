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

inline size_t repair_slot(const std::vector<Meta> &metas, const std::vector<Plan> &plan,
                          const Chain &chain, unsigned ranks, unsigned count,
                          unsigned anchor, uint64_t begin, uint64_t stride_bits) {
    if (!ranks || !count || anchor >= ranks || metas.size() != size_t(ranks) * count ||
        plan.size() != metas.size()) throw std::runtime_error("invalid repair job grid");
    if (chain.final || chain.cursor < begin || !stride_bits) return metas.size();
    const uint64_t slot = (chain.cursor - begin) / stride_bits;
    if (slot >= uint64_t(ranks) * count) return metas.size();
    const size_t id = size_t((anchor + slot / count) % ranks) * count + slot % count;
    const auto &m = metas[id];
    if (plan[id][0] || !m[PRESENT] || chain.cursor < m[BASE] || chain.cursor >= m[LIMIT])
        return metas.size();
    return id;
}

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

struct HistorySummary {
    uint64_t bytes;
    std::array<uint16_t, HISTORY> symbols;
};
static_assert(sizeof(HistorySummary) == sizeof(uint64_t) + HISTORY * sizeof(uint16_t),
              "history summary wire layout must be contiguous");

// Compose a job's tail with a symbolic incoming 32-KiB dictionary. References
// remain relative to the rank's incoming history, not to intermediate jobs.
inline void advance_symbols(const uint16_t *previous, const uint16_t *symbols,
                            size_t size, uint16_t *next) {
    const size_t n = std::min(size, HISTORY), keep = HISTORY - n;
    std::copy(previous + n, previous + HISTORY, next);
    for (size_t i = 0; i < n; ++i) {
        const unsigned symbol = symbols[size - n + i];
        if (symbol >= 256 && symbol < 65536U - HISTORY)
            throw std::runtime_error("invalid symbolic history marker");
        next[keep + i] = symbol < 256 ? uint16_t(symbol) : previous[symbol - (65536U - HISTORY)];
    }
}

struct SymbolSpan {
    const uint16_t *symbols;
    size_t size;
};

constexpr size_t HISTORY_TRACE_SLOTS = 65536;
struct HistoryTraceEntry {
    uint64_t position = UINT64_MAX;
    uint16_t symbol = 0;
};
inline size_t history_trace_slot(uint64_t position) {
    return size_t((position * UINT64_C(11400714819323198485)) >> 48);
}

inline void summarize_history(const std::vector<SymbolSpan> &spans, HistorySummary &summary) {
    std::vector<uint64_t> starts;
    std::vector<const uint16_t *> buffers;
    starts.reserve(spans.size()); buffers.reserve(spans.size());
    summary.bytes = 0;
    for (const auto &span : spans) {
        if (!span.size) continue;
        starts.push_back(summary.bytes); buffers.push_back(span.symbols);
        summary.bytes = add(summary.bytes, span.size);
    }
    const size_t count = size_t(std::min<uint64_t>(HISTORY, summary.bytes));
    const size_t keep = HISTORY - count;
    for (size_t i = 0; i < keep; ++i)
        summary.symbols[i] = uint16_t(65536U - HISTORY + count + i);
    if (!count) return;
    std::vector<HistoryTraceEntry> cache(HISTORY_TRACE_SLOTS);
    std::vector<uint64_t> path;
    path.reserve(starts.size());
    // DISCOVER has already expanded within-job references. Only markers in
    // the final window need tracing. Repeated quality bytes often share long
    // reference chains; cache resolved positions without changing job symbols.
    for (size_t i = 0; i < count; ++i) {
        uint64_t position = summary.bytes - count + i;
        size_t j = size_t(std::upper_bound(starts.begin(), starts.end(), position)
                          - starts.begin() - 1);
        uint16_t result;
        path.clear();
        for (;;) {
            const auto &entry = cache[history_trace_slot(position)];
            if (entry.position == position) {
                result = entry.symbol;
                break;
            }
            path.push_back(position);
            const unsigned symbol = buffers[j][size_t(position - starts[j])];
            if (symbol < 256) {
                result = uint16_t(symbol);
                break;
            }
            if (symbol < 65536U - HISTORY)
                throw std::runtime_error("invalid traced history marker");
            const unsigned distance = 65536U - symbol;
            if (distance > starts[j]) {
                result = uint16_t(65536U - (distance - starts[j]));
                break;
            }
            position = starts[j] - distance; // Strictly before the current job.
            do { --j; } while (starts[j] > position);
        }
        summary.symbols[keep + i] = result;
        // Direct-mapped replacement bounds memory; collisions only cause a
        // later cache miss, never reuse a result for a different position.
        for (uint64_t visited : path) {
            auto &entry = cache[history_trace_slot(visited)];
            entry.position = visited;
            entry.symbol = result;
        }
    }
}

inline Window resolve_history(const Window &previous, const HistorySummary &summary) {
    if (previous.size() > HISTORY) throw std::runtime_error("oversized summary dictionary");
    const size_t size = std::min<uint64_t>(HISTORY, previous.size() + std::min<uint64_t>(HISTORY, summary.bytes));
    Window next;
    next.reserve(size);
    for (size_t i = HISTORY - size; i < HISTORY; ++i) {
        unsigned symbol = summary.symbols[i];
        if (symbol >= 256) {
            if (symbol < 65536U - previous.size())
                throw std::runtime_error("invalid summary marker / missing history");
            symbol = previous[symbol - (65536U - previous.size())];
        }
        next.push_back(static_cast<unsigned char>(symbol));
    }
    return next;
}
} // namespace swbwa_mpi_gzip
#endif
