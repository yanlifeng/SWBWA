#include "fastq_frame.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace swbwa_fastq {
namespace {
Record validate(const unsigned char *data, size_t begin, const size_t *e, size_t cap)
{
    if (e[3] - begin > cap) throw std::runtime_error("FASTQ record exceeds cap");
    if (e[0] - begin < 3 || data[begin] != '@')
        throw std::runtime_error("invalid FASTQ header");
    if (e[2] - e[1] < 2 || data[e[1]] != '+')
        throw std::runtime_error("invalid FASTQ plus line");
    if (e[1] - e[0] != e[3] - e[2])
        throw std::runtime_error("FASTQ sequence/quality length mismatch");
    size_t name_end = begin + 1;
    while (name_end < e[0] - 1 && data[name_end] != ' ' && data[name_end] != '\t')
        ++name_end;
    if (name_end == begin + 1) throw std::runtime_error("empty FASTQ name token");
    return Record{begin, e[3], name_end - begin - 1};
}
void byte_check(unsigned char c)
{
    if (!c || c == '\r') throw std::runtime_error("FASTQ requires LF-only text without NUL/CR");
}
}

void LineOffsets::scan(const unsigned char *data, size_t size)
{
    if (size && !data) throw std::invalid_argument("null FASTQ data");
    ends.clear(); bytes = size;
    if (!size) return;
    if (std::memchr(data, 0, size) || std::memchr(data, '\r', size))
        throw std::runtime_error("FASTQ requires LF-only text without NUL/CR");
    const unsigned char *p = data, *end = data + size;
    while (p < end) {
        p = static_cast<const unsigned char *>(std::memchr(p, '\n', size_t(end - p)));
        if (!p) break;
        ends.push_back(size_t(++p - data));
    }
}

Split split(const LineOffsets &index, uint64_t before, bool after_lf)
{
    if (index.ends.size() > UINT64_MAX - before) throw std::overflow_error("FASTQ LF offset overflow");
    const bool boundary = !(before % 4) && after_lf;
    Split s = {0, 0, 0, 0, false, before / 4};
    const size_t first = 3 - size_t(before % 4);
    if (first >= index.ends.size()) {
        if (!boundary) s.prefix_end = s.body_begin = s.body_end = s.suffix_begin = index.bytes;
        return s;
    }
    const size_t last = first + ((index.ends.size() - 1 - first) / 4) * 4;
    if (!boundary) {
        s.prefix_end = s.body_begin = index.ends[first];
        s.prefix_completes = true;
    }
    s.body_end = s.suffix_begin = index.ends[last];
    return s;
}

Split split(const unsigned char *data, size_t size, uint64_t before,
            bool after_lf, LineOffsets &index)
{
    index.scan(data, size);
    return split(index, before, after_lf);
}

std::vector<Record> parse(const unsigned char *data, const LineOffsets &index,
                          size_t begin, size_t end, size_t cap)
{
    if (begin > end || end > index.bytes || !cap || (end && !data))
        throw std::invalid_argument("invalid FASTQ parse interval");
    std::vector<Record> records;
    auto first = std::upper_bound(index.ends.begin(), index.ends.end(), begin);
    auto last = std::upper_bound(index.ends.begin(), index.ends.end(), end);
    if ((last - first) % 4 || (begin != end && (last == first || *(last - 1) != end)))
        throw std::runtime_error("incomplete FASTQ record interval");
    records.reserve(size_t(last - first) / 4);
    while (first != last) {
        Record r = validate(data, begin, &*first, cap);
        records.push_back(r); begin = r.end; first += 4;
    }
    return records;
}

bool same_name(const unsigned char *a, const Record &x, const unsigned char *b, const Record &y)
{
    return x.name_size == y.name_size &&
        !std::memcmp(a + x.begin + 1, b + y.begin + 1, x.name_size);
}

OwnedChunk::OwnedChunk() : data(NULL), size(0), first_record(0), uncompressed_offset(0) {}
OwnedChunk::~OwnedChunk() { std::free(data); }
OwnedChunk::OwnedChunk(OwnedChunk &&other) noexcept : OwnedChunk() { *this = std::move(other); }
OwnedChunk &OwnedChunk::operator=(OwnedChunk &&other) noexcept
{
    if (this != &other) {
        std::free(data); data = other.release(); size = other.size;
        first_record = other.first_record; uncompressed_offset = other.uncompressed_offset;
        records = std::move(other.records);
    }
    return *this;
}
char *OwnedChunk::release() { char *p = data; data = NULL; return p; }

Frame::Frame(size_t cap, uint64_t offset, uint64_t first) :
    cap_(cap), head_(0), taken_(0), scan_(0), record_begin_(0), line_begin_(0), phase_(0),
    base_(offset), ordinal_(first), failed_(false), finished_(false)
{
    if (!cap) throw std::invalid_argument("zero FASTQ record cap");
}
size_t Frame::available() const { return records_.size() - taken_; }
size_t Frame::partial_bytes() const { return data_.size() - record_begin_; }

void Frame::compact()
{
    if (head_ && head_ >= data_.size() / 2) {
        data_.erase(data_.begin(), data_.begin() + head_);
        base_ += head_; scan_ -= head_; record_begin_ -= head_; line_begin_ -= head_;
        for (size_t i = 0; i < phase_; ++i) line_ends_[i] -= head_;
        for (size_t i = taken_; i < records_.size(); ++i) {
            records_[i].begin -= head_; records_[i].end -= head_;
        }
        head_ = 0;
    }
    if (taken_ && taken_ >= records_.size() / 2) {
        records_.erase(records_.begin(), records_.begin() + taken_); taken_ = 0;
    }
}

void Frame::append(const void *input, size_t size)
{
    if (failed_ || finished_) throw std::logic_error("FASTQ frame not appendable");
    try {
        if (size && !input) throw std::invalid_argument("null FASTQ append");
        compact();
        if (size > UINT64_MAX - base_ - data_.size() || size > data_.max_size() - data_.size())
            throw std::overflow_error("FASTQ byte offset overflow");
        const unsigned char *p = static_cast<const unsigned char *>(input);
        if (size) data_.insert(data_.end(), p, p + size);
        for (; scan_ < data_.size(); ++scan_) {
            byte_check(data_[scan_]);
            if (scan_ + 1 - record_begin_ > cap_) throw std::runtime_error("FASTQ record exceeds cap");
            if (data_[scan_] != '\n') continue;
            line_ends_[phase_++] = scan_ + 1; line_begin_ = scan_ + 1;
            if (phase_ == 4) {
                records_.push_back(validate(data_.data(), record_begin_, line_ends_, cap_));
                record_begin_ = scan_ + 1; phase_ = 0;
            }
        }
    } catch (...) { failed_ = true; throw; }
}

void Frame::finish()
{
    if (failed_) throw std::logic_error("FASTQ frame failed");
    if (partial_bytes()) {
        failed_ = true;
        throw std::runtime_error("truncated FASTQ: four LF-terminated lines required");
    }
    finished_ = true;
}

OwnedChunk Frame::take_prefix(size_t count)
{
    if (failed_) throw std::logic_error("FASTQ frame failed");
    if (count > available()) throw std::out_of_range("FASTQ prefix exceeds available records");
    if (count > UINT64_MAX - ordinal_) throw std::overflow_error("FASTQ record ordinal overflow");
    OwnedChunk out;
    out.first_record = ordinal_; out.uncompressed_offset = base_ + head_;
    if (!count) return out;
    const size_t end = records_[taken_ + count - 1].end;
    out.size = end - head_;
    out.data = static_cast<char *>(std::malloc(out.size + 1));
    if (!out.data) throw std::bad_alloc();
    std::memcpy(out.data, data_.data() + head_, out.size); out.data[out.size] = 0;
    out.records.reserve(count);
    for (size_t i = taken_; i < taken_ + count; ++i) {
        Record r = records_[i]; r.begin -= head_; r.end -= head_; out.records.push_back(r);
    }
    head_ = end; taken_ += count; ordinal_ += count;
    compact();
    return out;
}
}
