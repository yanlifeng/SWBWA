#include "swbwa_gzip.h"
#include "gzip_stream.h"
#include "fastq_layout.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
thread_local char last_error[512] = "";

double now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

void set_error(char *destination, const char *message) {
    std::snprintf(destination, 512, "%s", message);
}

struct Free {
    void operator()(char *p) const { std::free(p); }
};
using Owned = std::unique_ptr<char, Free>;

struct Segment {
    SwbwaDecodedView view;
    size_t offset = 0, consumed_lines = 0;
};

struct Input {
    SwbwaGzipReader gzip;
    std::deque<Segment> pending;
    size_t batch = 16U << 20, bytes = 0, lines = 0;
    bool eof = false;
    swbwa_gzip_input_metrics metrics{};
    SwbwaLineIndex delivered_index;

    Input(const char *path, size_t size, unsigned workers, const char *backend)
        : gzip(path, size, workers, backend, nullptr, true) {}

    void prepare(size_t target_records = 0) {
        while ((target_records ? lines / 4 < target_records : (bytes < batch || lines < 4)) && !eof) {
            Segment segment;
            auto &v = segment.view;
            ++metrics.refills;
            const size_t request = target_records ? size_t(1U << 20) :
                (bytes < batch ? batch - bytes : size_t(64U << 10));
            if (gzip.supports_views()) {
                if (gzip.peek(request, v)) gzip.consume(v.size);
            } else {
                const size_t capacity = std::min(request, size_t(1U << 20));
                v.owner = std::shared_ptr<unsigned char>(new unsigned char[capacity],
                    std::default_delete<unsigned char[]>());
                v.data = v.owner.get();
                v.size = gzip.read(v.data, capacity, &v.lines);
            }
            if (!v.size) {
                eof = true;
                if (lines % 4 || (bytes && pending.back().view.data[pending.back().view.size - 1] != '\n'))
                    throw std::runtime_error("truncated FASTQ: expected four LF-terminated lines (final LF required)");
                break;
            }
            if (v.size > size_t(std::numeric_limits<int64_t>::max()) - bytes)
                throw std::length_error("FASTQ record exceeds int64 byte range");
            bytes += v.size;
            lines += v.lines.lines();
            pending.push_back(std::move(segment));
            // rapid may know EOF before the caller drains its last view.
            if (gzip.eof()) {
                eof = true;
                if (lines % 4 || pending.back().view.data[pending.back().view.size - 1] != '\n')
                    throw std::runtime_error("truncated FASTQ: expected four LF-terminated lines (final LF required)");
            }
        }
    }

    size_t prefix_bytes(size_t count) const {
        size_t total = 0;
        for (const auto &segment : pending) {
            const auto &v = segment.view;
            const size_t available = v.lines.lines() - segment.consumed_lines;
            if (count <= available) {
                return total + v.lines.end_of_line(v.data, segment.consumed_lines + count) - segment.offset;
            }
            total += v.size - segment.offset;
            count -= available;
        }
        throw std::logic_error("FASTQ prefix exceeds buffered line count");
    }

    size_t ready_records() const {
        if (!lines) return 0;
        size_t remaining = batch, count = 0;
        for (const auto &segment : pending) {
            const auto &v = segment.view;
            const size_t available = v.size - segment.offset;
            if (remaining >= available) {
                count += v.lines.lines() - segment.consumed_lines;
                remaining -= available;
            } else {
                count += v.lines.lines_before(v.data, segment.offset + remaining) - segment.consumed_lines;
                break;
            }
        }
        return std::max<size_t>(1, count / 4);
    }

    Owned extract(size_t records, size_t &size) {
        size = prefix_bytes(records * 4);
        if (size == std::numeric_limits<size_t>::max()) throw std::bad_alloc();
        Owned result(static_cast<char *>(std::malloc(size + 1)));
        if (!result) throw std::bad_alloc();
        delivered_index.clear();
        size_t copied = 0;
        const double start = now();
        while (copied < size) {
            auto &segment = pending.front();
            const auto &v = segment.view;
            const size_t n = std::min(size - copied, v.size - segment.offset);
            std::memcpy(result.get() + copied, v.data + segment.offset, n);
            const size_t previous = delivered_index.lines();
            delivered_index.append_slice(v.lines, v.data + segment.offset, segment.offset, n);
            segment.consumed_lines += delivered_index.lines() - previous;
            segment.offset += n;
            copied += n;
            ++metrics.copy_calls;
            if (segment.offset == v.size) pending.pop_front();
        }
        metrics.copy_s += now() - start;
        metrics.copy_bytes += copied;
        bytes -= size;
        lines -= records * 4;
        result.get()[size] = '\0';
        return result;
    }

    swbwa_gzip_input_metrics snapshot() const {
        auto result = metrics;
        const auto &m = gzip.metrics();
        result.input_bytes = m.input_bytes;
        result.decoded_bytes = m.output_bytes;
        result.jobs = m.jobs;
        result.accepted_jobs = m.accepted_jobs;
        result.fallback_bytes = m.fallback_bytes;
        result.capacity_retries = m.capacity_retries;
        result.source_read_bytes = m.source_read_bytes;
        result.prefetch_bytes = m.prefetch_bytes;
        result.peak_symbol_bytes = m.peak_symbol_bytes;
        result.read_s = m.read_s;
        result.kernel_s = m.kernel_s;
        result.queue_s = m.queue_s;
        result.serial_s = m.serial_s;
        result.line_index_s = m.line_index_s;
        result.io_s = m.io_s;
        result.discover_run_s = m.discover_run_s;
        result.resolve_run_s = m.resolve_run_s;
        result.repair_s = m.repair_s;
        result.input_prepare_s = m.input_prepare_s;
        result.task_build_s = m.task_build_s;
        result.prefetch_io_s = m.prefetch_io_s;
        result.discover_wait_s = m.discover_wait_s;
        result.resolve_wait_s = m.resolve_wait_s;
        return result;
    }
};

}

struct swbwa_gzip_reader {
    std::unique_ptr<Input> input[2];
    char error[512] = "";
    uint64_t chunks = 0, records = 0, last_records = 0;
    bool started = false;
};

extern "C" swbwa_gzip_reader *swbwa_gzip_open(const char *path1, const char *path2,
                                            size_t chunk, unsigned workers,
                                            const char *backend) {
    last_error[0] = '\0';
    try {
        if (!path1 || !*path1 || (path2 && !*path2) || !backend ||
            (std::strcmp(backend, "rapid") && std::strcmp(backend, "zlib")))
            throw std::invalid_argument("gzip open requires paths and backend rapid or zlib");
        std::unique_ptr<swbwa_gzip_reader> reader(new swbwa_gzip_reader);
        reader->input[0].reset(new Input(path1, chunk, workers, backend));
        if (path2) reader->input[1].reset(new Input(path2, chunk, workers, backend));
        return reader.release();
    } catch (const std::exception &e) { set_error(last_error, e.what()); }
    catch (...) { set_error(last_error, "unknown gzip open failure"); }
    return nullptr;
}

extern "C" int swbwa_gzip_set_batch_bytes(swbwa_gzip_reader *reader, size_t bytes) {
    if (!reader) { set_error(last_error, "null gzip reader"); return -1; }
    if (reader->error[0]) return -1;
    if (reader->started || !bytes || bytes >= (uint64_t)std::numeric_limits<int64_t>::max()) {
        set_error(reader->error, "set_batch_bytes requires positive size below INT64_MAX before first read");
        return -1;
    }
    for (auto &input : reader->input) if (input) input->batch = bytes;
    return 0;
}

extern "C" int swbwa_gzip_read_next(swbwa_gzip_reader *reader,
                                   char **r1, int64_t *n1, char **r2, int64_t *n2) {
    if (r1) *r1 = nullptr;
    if (r2) *r2 = nullptr;
    if (n1) *n1 = 0;
    if (n2) *n2 = 0;
    if (reader) reader->last_records = 0;
    if (!reader || !r1 || !r2 || !n1 || !n2 || r1 == r2 || n1 == n2) {
        set_error(reader ? reader->error : last_error, "invalid gzip read_next arguments");
        return -1;
    }
    if (reader->error[0]) return -1;
    reader->started = true;
    try {
        auto &first = *reader->input[0];
        first.prepare();
        size_t records = first.ready_records();
        if (reader->input[1]) {
            auto &second = *reader->input[1];
            second.prepare(records ? records : 1);
            if ((!records && second.lines) || second.lines / 4 < records)
                throw std::runtime_error("paired FASTQ record count mismatch");
        }
        if (!records) return 0;
        size_t size1 = 0, size2 = 0;
        Owned one = first.extract(records, size1), two;
        if (reader->input[1]) two = reader->input[1]->extract(records, size2);
        *n1 = (int64_t)size1; *n2 = (int64_t)size2;
        *r1 = one.release(); *r2 = two.release();
        ++reader->chunks;
        reader->records += records;
        reader->last_records = records;
        first.metrics.delivered_bytes += size1;
        if (reader->input[1]) reader->input[1]->metrics.delivered_bytes += size2;
        return 1;
    } catch (const std::exception &e) { set_error(reader->error, e.what()); }
    catch (...) { set_error(reader->error, "unknown gzip read failure"); }
    return -1;
}

extern "C" int swbwa_gzip_next(swbwa_gzip_reader *reader,
                              char **r1, int64_t *n1, char **r2, int64_t *n2) {
    return swbwa_gzip_read_next(reader, r1, n1, r2, n2);
}

extern "C" const char *swbwa_gzip_error(const swbwa_gzip_reader *reader) {
    return reader ? reader->error : last_error;
}

extern "C" void swbwa_gzip_close(swbwa_gzip_reader *reader) { delete reader; }

extern "C" int swbwa_gzip_check_pe_layout(swbwa_gzip_reader *reader,
                                         const char *r1, int64_t n1,
                                         const char *r2, int64_t n2,
                                         unsigned count) {
    if (!reader) { set_error(last_error, "null gzip reader"); return -1; }
    if (reader->error[0]) return -1;
    try {
        if (!count || count > 384) throw std::invalid_argument("CPE count must be 1..384");
        if (!reader->input[1]) return 0;
        const auto &a = reader->input[0]->delivered_index;
        const auto &b = reader->input[1]->delivered_index;
        if (!reader->last_records || !r1 || !r2 || n1 <= 0 || n2 <= 0 ||
            (uint64_t)n1 != a.bytes() || (uint64_t)n2 != b.bytes())
            throw std::invalid_argument("layout check requires the latest unmodified PE batch");
        const auto *data1 = reinterpret_cast<const unsigned char *>(r1);
        const auto *data2 = reinterpret_cast<const unsigned char *>(r2);
        swbwa_fastq_check_layout(data1, a, data2, b, count);
        return 0;
    } catch (const std::exception &e) { set_error(reader->error, e.what()); }
    catch (...) { set_error(reader->error, "unknown PE layout failure"); }
    return -1;
}

extern "C" int swbwa_gzip_get_metrics(const swbwa_gzip_reader *reader,
                                     swbwa_gzip_metrics *metrics) {
    if (!reader || !metrics) { set_error(last_error, "invalid gzip metrics arguments"); return -1; }
    *metrics = swbwa_gzip_metrics{};
    metrics->chunks = reader->chunks;
    metrics->records = reader->records;
    metrics->last_records = reader->last_records;
    metrics->paired = reader->input[1] != nullptr;
    for (int i = 0; i < 2; ++i)
        if (reader->input[i]) metrics->input[i] = reader->input[i]->snapshot();
    return 0;
}

extern "C" int swbwa_gzip_report(const swbwa_gzip_reader *reader, FILE *stream) {
    swbwa_gzip_metrics metrics;
    if (!stream || swbwa_gzip_get_metrics(reader, &metrics)) {
        set_error(last_error, "invalid gzip report arguments"); return -1;
    }
    if (std::fprintf(stream, "[SWBWA gzip] paired=%d chunks=%llu records=%llu\n",
                     metrics.paired, (unsigned long long)metrics.chunks,
                     (unsigned long long)metrics.records) < 0) return -1;
    for (int i = 0; i <= metrics.paired; ++i) {
        const auto &m = metrics.input[i];
        if (std::fprintf(stream,
            "[SWBWA gzip input] mate=%d compressed=%llu decoded=%llu delivered=%llu "
            "jobs=%llu accepted=%llu fallback=%llu copy_bytes=%llu "
            "read_s=%.6f kernel_s=%.6f queue_s=%.6f serial_s=%.6f index_s=%.6f copy_s=%.6f "
            "io_s=%.6f discover_run_s=%.6f resolve_run_s=%.6f repair_s=%.6f "
            "input_prepare_s=%.6f task_build_s=%.6f prefetch_io_s=%.6f "
            "discover_wait_s=%.6f resolve_wait_s=%.6f\n",
            i + 1, (unsigned long long)m.input_bytes, (unsigned long long)m.decoded_bytes,
            (unsigned long long)m.delivered_bytes, (unsigned long long)m.jobs,
            (unsigned long long)m.accepted_jobs, (unsigned long long)m.fallback_bytes,
            (unsigned long long)m.copy_bytes, m.read_s, m.kernel_s, m.queue_s,
            m.serial_s, m.line_index_s, m.copy_s, m.io_s, m.discover_run_s,
            m.resolve_run_s, m.repair_s, m.input_prepare_s, m.task_build_s,
            m.prefetch_io_s, m.discover_wait_s, m.resolve_wait_s) < 0) return -1;
    }
    return 0;
}
