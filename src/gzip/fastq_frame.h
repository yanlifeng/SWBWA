#ifndef SWBWA_FASTQ_FRAME_H
#define SWBWA_FASTQ_FRAME_H

#include <cstddef>
#include <cstdint>
#include <vector>

namespace swbwa_fastq {

static const size_t RECORD_CAP = 1U << 20;

struct Record {
    size_t begin, end, name_size; // Name starts at begin + 1; exact token.
};

struct LineOffsets {
    size_t bytes;
    std::vector<size_t> ends; // Exclusive LF offsets relative to this extent.
    LineOffsets() : bytes(0) {}
    void scan(const unsigned char *data, size_t size);
};

struct Split {
    size_t prefix_end, body_begin, body_end, suffix_begin;
    bool prefix_completes;
    uint64_t first_record; // First record ending here, if any.
};

/* global_lf_before counts LFs strictly before data[0]. starts_after_lf is true
 * at file origin or when its preceding global byte is LF. No '@' search.
 * A mid-record extent with no fourth LF is prefix only, not also suffix.
 */
Split split(const LineOffsets &, uint64_t global_lf_before, bool starts_after_lf);
Split split(const unsigned char *, size_t, uint64_t global_lf_before,
            bool starts_after_lf, LineOffsets &);

std::vector<Record> parse(const unsigned char *, const LineOffsets &,
                          size_t begin, size_t end, size_t cap = RECORD_CAP);
bool same_name(const unsigned char *, const Record &,
               const unsigned char *, const Record &);

class OwnedChunk {
public:
    char *data;
    size_t size;
    uint64_t first_record, uncompressed_offset;
    std::vector<Record> records;
    OwnedChunk();
    ~OwnedChunk();
    OwnedChunk(OwnedChunk &&) noexcept;
    OwnedChunk &operator=(OwnedChunk &&) noexcept;
    char *release(); // Caller subsequently free()s it; metadata remains valid.
private:
    OwnedChunk(const OwnedChunk &) = delete;
    OwnedChunk &operator=(const OwnedChunk &) = delete;
};

/* Ordered incremental stream. Scans each incoming byte once; bounded partial
 * record, block allocations rather than per-record allocations. Complete
 * records stay buffered until take_prefix. Exceptions make append/finish
 * sticky-failed. Caller bounds the number of complete records it queues.
 */
class Frame {
public:
    explicit Frame(size_t cap = RECORD_CAP, uint64_t byte_offset = 0,
                   uint64_t first_record = 0);
    void append(const void *data, size_t size);
    void finish(); // Requires final LF and no partial record.
    size_t available() const;
    size_t partial_bytes() const;
    OwnedChunk take_prefix(size_t records);
private:
    void compact();
    size_t cap_, head_, taken_, scan_, record_begin_, line_begin_, phase_;
    size_t line_ends_[4];
    uint64_t base_, ordinal_;
    bool failed_, finished_;
    std::vector<unsigned char> data_;
    std::vector<Record> records_;
};

}
#endif
