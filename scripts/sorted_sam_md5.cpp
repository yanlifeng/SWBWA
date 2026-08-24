#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct SamRecord {
    std::uint64_t id;
    std::string line;
};

class Md5 {
public:
    Md5()
        : state_{{0x67452301U, 0xefcdab89U, 0x98badcfeU, 0x10325476U}},
          total_size_(0), buffer_size_(0)
    {
    }

    void update(const void *data, std::size_t size)
    {
        const std::uint8_t *bytes =
            static_cast<const std::uint8_t *>(data);

        total_size_ += size;
        if (buffer_size_ != 0) {
            const std::size_t available = buffer_.size() - buffer_size_;
            const std::size_t copied = size < available ? size : available;

            std::copy(bytes, bytes + copied, buffer_.begin() + buffer_size_);
            buffer_size_ += copied;
            bytes += copied;
            size -= copied;
            if (buffer_size_ != buffer_.size()) return;
            transform(buffer_.data());
            buffer_size_ = 0;
        }
        while (size >= buffer_.size()) {
            transform(bytes);
            bytes += buffer_.size();
            size -= buffer_.size();
        }
        if (size != 0) {
            std::copy(bytes, bytes + size, buffer_.begin());
            buffer_size_ = size;
        }
    }

    std::array<std::uint8_t, 16> finish()
    {
        const std::uint64_t bit_size = total_size_ * 8;
        const std::uint8_t padding[64] = {0x80};
        std::uint8_t encoded_size[8];
        const std::size_t padding_size = buffer_size_ < 56
                                       ? 56 - buffer_size_
                                       : 120 - buffer_size_;

        for (unsigned i = 0; i < 8; ++i)
            encoded_size[i] =
                static_cast<std::uint8_t>(bit_size >> (i * 8));
        update(padding, padding_size);
        update(encoded_size, sizeof(encoded_size));

        std::array<std::uint8_t, 16> digest;
        for (unsigned i = 0; i < state_.size(); ++i) {
            for (unsigned byte = 0; byte < 4; ++byte) {
                digest[i * 4 + byte] = static_cast<std::uint8_t>(
                    state_[i] >> (byte * 8));
            }
        }
        return digest;
    }

private:
    static std::uint32_t rotate_left(std::uint32_t value, unsigned shift)
    {
        return (value << shift) | (value >> (32 - shift));
    }

    void transform(const std::uint8_t block[64])
    {
        std::uint32_t words[16];
        for (unsigned i = 0; i < 16; ++i) {
            words[i] = static_cast<std::uint32_t>(block[i * 4]) |
                       (static_cast<std::uint32_t>(block[i * 4 + 1]) << 8) |
                       (static_cast<std::uint32_t>(block[i * 4 + 2]) << 16) |
                       (static_cast<std::uint32_t>(block[i * 4 + 3]) << 24);
        }

        std::uint32_t a = state_[0];
        std::uint32_t b = state_[1];
        std::uint32_t c = state_[2];
        std::uint32_t d = state_[3];

#define MD5_F(x, y, z) (((x) & (y)) | (~(x) & (z)))
#define MD5_G(x, y, z) (((x) & (z)) | ((y) & ~(z)))
#define MD5_H(x, y, z) ((x) ^ (y) ^ (z))
#define MD5_I(x, y, z) ((y) ^ ((x) | ~(z)))
#define MD5_STEP(function, aa, bb, cc, dd, word, constant, shift) do { \
    (aa) += function((bb), (cc), (dd)) + (word) + (constant); \
    (aa) = rotate_left((aa), (shift)) + (bb); \
} while (0)

        MD5_STEP(MD5_F, a, b, c, d, words[0],  0xd76aa478U, 7);
        MD5_STEP(MD5_F, d, a, b, c, words[1],  0xe8c7b756U, 12);
        MD5_STEP(MD5_F, c, d, a, b, words[2],  0x242070dbU, 17);
        MD5_STEP(MD5_F, b, c, d, a, words[3],  0xc1bdceeeU, 22);
        MD5_STEP(MD5_F, a, b, c, d, words[4],  0xf57c0fafU, 7);
        MD5_STEP(MD5_F, d, a, b, c, words[5],  0x4787c62aU, 12);
        MD5_STEP(MD5_F, c, d, a, b, words[6],  0xa8304613U, 17);
        MD5_STEP(MD5_F, b, c, d, a, words[7],  0xfd469501U, 22);
        MD5_STEP(MD5_F, a, b, c, d, words[8],  0x698098d8U, 7);
        MD5_STEP(MD5_F, d, a, b, c, words[9],  0x8b44f7afU, 12);
        MD5_STEP(MD5_F, c, d, a, b, words[10], 0xffff5bb1U, 17);
        MD5_STEP(MD5_F, b, c, d, a, words[11], 0x895cd7beU, 22);
        MD5_STEP(MD5_F, a, b, c, d, words[12], 0x6b901122U, 7);
        MD5_STEP(MD5_F, d, a, b, c, words[13], 0xfd987193U, 12);
        MD5_STEP(MD5_F, c, d, a, b, words[14], 0xa679438eU, 17);
        MD5_STEP(MD5_F, b, c, d, a, words[15], 0x49b40821U, 22);

        MD5_STEP(MD5_G, a, b, c, d, words[1],  0xf61e2562U, 5);
        MD5_STEP(MD5_G, d, a, b, c, words[6],  0xc040b340U, 9);
        MD5_STEP(MD5_G, c, d, a, b, words[11], 0x265e5a51U, 14);
        MD5_STEP(MD5_G, b, c, d, a, words[0],  0xe9b6c7aaU, 20);
        MD5_STEP(MD5_G, a, b, c, d, words[5],  0xd62f105dU, 5);
        MD5_STEP(MD5_G, d, a, b, c, words[10], 0x02441453U, 9);
        MD5_STEP(MD5_G, c, d, a, b, words[15], 0xd8a1e681U, 14);
        MD5_STEP(MD5_G, b, c, d, a, words[4],  0xe7d3fbc8U, 20);
        MD5_STEP(MD5_G, a, b, c, d, words[9],  0x21e1cde6U, 5);
        MD5_STEP(MD5_G, d, a, b, c, words[14], 0xc33707d6U, 9);
        MD5_STEP(MD5_G, c, d, a, b, words[3],  0xf4d50d87U, 14);
        MD5_STEP(MD5_G, b, c, d, a, words[8],  0x455a14edU, 20);
        MD5_STEP(MD5_G, a, b, c, d, words[13], 0xa9e3e905U, 5);
        MD5_STEP(MD5_G, d, a, b, c, words[2],  0xfcefa3f8U, 9);
        MD5_STEP(MD5_G, c, d, a, b, words[7],  0x676f02d9U, 14);
        MD5_STEP(MD5_G, b, c, d, a, words[12], 0x8d2a4c8aU, 20);

        MD5_STEP(MD5_H, a, b, c, d, words[5],  0xfffa3942U, 4);
        MD5_STEP(MD5_H, d, a, b, c, words[8],  0x8771f681U, 11);
        MD5_STEP(MD5_H, c, d, a, b, words[11], 0x6d9d6122U, 16);
        MD5_STEP(MD5_H, b, c, d, a, words[14], 0xfde5380cU, 23);
        MD5_STEP(MD5_H, a, b, c, d, words[1],  0xa4beea44U, 4);
        MD5_STEP(MD5_H, d, a, b, c, words[4],  0x4bdecfa9U, 11);
        MD5_STEP(MD5_H, c, d, a, b, words[7],  0xf6bb4b60U, 16);
        MD5_STEP(MD5_H, b, c, d, a, words[10], 0xbebfbc70U, 23);
        MD5_STEP(MD5_H, a, b, c, d, words[13], 0x289b7ec6U, 4);
        MD5_STEP(MD5_H, d, a, b, c, words[0],  0xeaa127faU, 11);
        MD5_STEP(MD5_H, c, d, a, b, words[3],  0xd4ef3085U, 16);
        MD5_STEP(MD5_H, b, c, d, a, words[6],  0x04881d05U, 23);
        MD5_STEP(MD5_H, a, b, c, d, words[9],  0xd9d4d039U, 4);
        MD5_STEP(MD5_H, d, a, b, c, words[12], 0xe6db99e5U, 11);
        MD5_STEP(MD5_H, c, d, a, b, words[15], 0x1fa27cf8U, 16);
        MD5_STEP(MD5_H, b, c, d, a, words[2],  0xc4ac5665U, 23);

        MD5_STEP(MD5_I, a, b, c, d, words[0],  0xf4292244U, 6);
        MD5_STEP(MD5_I, d, a, b, c, words[7],  0x432aff97U, 10);
        MD5_STEP(MD5_I, c, d, a, b, words[14], 0xab9423a7U, 15);
        MD5_STEP(MD5_I, b, c, d, a, words[5],  0xfc93a039U, 21);
        MD5_STEP(MD5_I, a, b, c, d, words[12], 0x655b59c3U, 6);
        MD5_STEP(MD5_I, d, a, b, c, words[3],  0x8f0ccc92U, 10);
        MD5_STEP(MD5_I, c, d, a, b, words[10], 0xffeff47dU, 15);
        MD5_STEP(MD5_I, b, c, d, a, words[1],  0x85845dd1U, 21);
        MD5_STEP(MD5_I, a, b, c, d, words[8],  0x6fa87e4fU, 6);
        MD5_STEP(MD5_I, d, a, b, c, words[15], 0xfe2ce6e0U, 10);
        MD5_STEP(MD5_I, c, d, a, b, words[6],  0xa3014314U, 15);
        MD5_STEP(MD5_I, b, c, d, a, words[13], 0x4e0811a1U, 21);
        MD5_STEP(MD5_I, a, b, c, d, words[4],  0xf7537e82U, 6);
        MD5_STEP(MD5_I, d, a, b, c, words[11], 0xbd3af235U, 10);
        MD5_STEP(MD5_I, c, d, a, b, words[2],  0x2ad7d2bbU, 15);
        MD5_STEP(MD5_I, b, c, d, a, words[9],  0xeb86d391U, 21);

#undef MD5_STEP
#undef MD5_I
#undef MD5_H
#undef MD5_G
#undef MD5_F
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
    }

    std::array<std::uint32_t, 4> state_;
    std::uint64_t total_size_;
    std::array<std::uint8_t, 64> buffer_;
    std::size_t buffer_size_;
};

static std::uint64_t parse_name_id(const std::string &line,
                                   std::size_t line_no)
{
    const std::size_t qname_end = line.find_first_of("\t ");
    const std::size_t end = qname_end == std::string::npos
                          ? line.size() : qname_end;

    if (end == 0)
        throw std::runtime_error("empty QNAME at line " +
                                 std::to_string(line_no));

    const std::size_t dot = line.rfind('.', end - 1);
    if (dot == std::string::npos || dot + 1 >= end)
        throw std::runtime_error("missing numeric QNAME suffix at line " +
                                 std::to_string(line_no));

    std::uint64_t value = 0;
    std::size_t pos = dot + 1;
    bool has_digit = false;
    while (pos < end && line[pos] >= '0' && line[pos] <= '9') {
        const unsigned digit = static_cast<unsigned>(line[pos] - '0');
        has_digit = true;
        if (value > (std::numeric_limits<std::uint64_t>::max() - digit) / 10)
            throw std::runtime_error("QNAME suffix overflow at line " +
                                     std::to_string(line_no));
        value = value * 10 + digit;
        ++pos;
    }
    if (!has_digit)
        throw std::runtime_error("empty numeric QNAME suffix at line " +
                                 std::to_string(line_no));
    return value;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " RESULT.md5 INPUT.sam [INPUT.sam ...]\n";
        return 1;
    }

    try {
        std::vector<SamRecord> records;
        std::size_t line_no = 0;
        const std::size_t initial_record_capacity =
            static_cast<std::size_t>(1) << 25;

        records.reserve(initial_record_capacity);
        std::cerr << "Reserved space for " << initial_record_capacity
                  << " SAM records\n";

        for (int arg = 2; arg < argc; ++arg) {
            std::vector<char> input_buffer(8U * 1024U * 1024U);
            std::ifstream input;
            input.rdbuf()->pubsetbuf(input_buffer.data(), input_buffer.size());
            input.open(argv[arg], std::ios::binary);
            if (!input)
                throw std::runtime_error("cannot open input: " +
                                         std::string(argv[arg]));

            std::string line;
            while (std::getline(input, line)) {
                ++line_no;
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                records.push_back(SamRecord{
                    parse_name_id(line, line_no), std::move(line)
                });
                if (line_no % 1000000 == 0)
                    std::cerr << "Loaded " << line_no << " SAM records\n";
            }
            if (!input.eof())
                throw std::runtime_error("error while reading input: " +
                                         std::string(argv[arg]));
        }

        std::cerr << "Finished loading " << records.size()
                  << " SAM records\nSorting records by ID...\n";
        std::stable_sort(records.begin(), records.end(),
            [](const SamRecord &a, const SamRecord &b) {
                return a.id < b.id;
            });
        std::cerr << "Sorting finished\nHashing sorted SAM stream...\n";

        Md5 md5;
        const char newline = '\n';

        for (std::size_t i = 0; i < records.size(); ++i) {
            const std::string &line = records[i].line;
            if (!line.empty())
                md5.update(line.data(), line.size());
            md5.update(&newline, 1);
        }

        const std::array<std::uint8_t, 16> digest = md5.finish();
        static const char hex[] = "0123456789abcdef";
        std::ofstream result(argv[1], std::ios::binary | std::ios::trunc);
        if (!result)
            throw std::runtime_error("cannot open result: " +
                                     std::string(argv[1]));
        for (std::size_t i = 0; i < digest.size(); ++i)
            result << hex[digest[i] >> 4] << hex[digest[i] & 15];
        result << '\n';
        if (!result)
            throw std::runtime_error("failed to write result: " +
                                     std::string(argv[1]));

        std::cout << "Hashed " << records.size()
                  << " sorted SAM records\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 2;
    }
}
