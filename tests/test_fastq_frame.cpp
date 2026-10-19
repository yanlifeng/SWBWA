#include "fastq_frame.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fq = swbwa_fastq;
template<class F> static void rejects(F f)
{
    bool caught = false;
    try { f(); } catch (const std::exception &) { caught = true; }
    assert(caught);
}
static const unsigned char *ptr(const std::string &s)
{
    return reinterpret_cast<const unsigned char *>(s.data());
}
static void stream(const std::string &s, size_t a, size_t b)
{
    fq::Frame f(fq::RECORD_CAP, 123, 7);
    f.append(s.data(), a); f.append(s.data() + a, b - a); f.append(s.data() + b, s.size() - b);
    f.finish(); assert(f.available() == 3 && !f.partial_bytes());
    std::string joined;
    for (size_t i = 0; i < 3; ++i) {
        auto c = f.take_prefix(1);
        assert(c.first_record == 7 + i && c.uncompressed_offset == 123 + joined.size());
        assert(c.data[c.size] == 0 && c.records.size() == 1);
        assert(c.records[0].begin == 0 && c.records[0].end == c.size);
        joined.append(c.data, c.size);
    }
    assert(joined == s && !f.available());
    rejects([&]() { f.take_prefix(1); });
    rejects([&]() { f.append("", 0); });
}
static void split_check(const std::string &s, const std::vector<size_t> &cuts)
{
    size_t start = 0, lfs = 0, ordinal = 0;
    std::string carry, joined;
    for (size_t end : cuts) {
        fq::LineOffsets index;
        const auto split = fq::split(ptr(s) + start, end - start, lfs,
                                    !start || s[start - 1] == '\n', index);
        if (split.prefix_end) carry.append(s, start, split.prefix_end);
        size_t ended = 0;
        if (split.prefix_completes) {
            fq::Frame one;
            one.append(carry.data(), carry.size()); one.finish();
            assert(one.available() == 1 && split.first_record == ordinal);
            joined += carry; carry.clear(); ++ordinal; ++ended;
        }
        const auto records = fq::parse(ptr(s) + start, index, split.body_begin, split.body_end);
        for (const auto &r : records) {
            assert(split.first_record + ended == ordinal);
            joined.append(s, start + r.begin, r.end - r.begin); ++ordinal; ++ended;
        }
        if (!split.prefix_end || split.prefix_completes)
            carry.append(s, start + split.suffix_begin, end - start - split.suffix_begin);
        assert(carry.size() <= fq::RECORD_CAP);
        lfs += index.ends.size(); start = end;
    }
    assert(start == s.size() && carry.empty() && ordinal == 3 && joined == s);
}
int main()
{
    const std::string text = "@one comment\nAC\n+text\n@+\n@two\tcomment\nT\n+\n+\n@three\n\n+\n\n";
    for (size_t a = 0; a <= text.size(); ++a)
        for (size_t b = a; b <= text.size(); ++b) {
            stream(text, a, b); split_check(text, {a, b, text.size()});
        }
    std::vector<size_t> bytes;
    fq::Frame f;
    std::string joined;
    for (size_t i = 0; i < text.size(); ++i) {
        bytes.push_back(i + 1); f.append(text.data() + i, 1);
        if (f.available()) { auto c = f.take_prefix(1); joined.append(c.data, c.size); }
    }
    f.finish(); assert(joined == text); split_check(text, bytes);
    for (size_t n = 0; n <= 3; ++n) {
        fq::Frame all; all.append(text.data(), text.size());
        auto a = all.take_prefix(n); auto b = all.take_prefix(3 - n);
        assert(std::string(a.data ? a.data : "", a.size) + std::string(b.data ? b.data : "", b.size) == text);
        all.finish();
    }
    const std::vector<std::string> bad = {
        "a\nA\n+\nI\n", "@\nA\n+\nI\n", "@ x\nA\n+\nI\n",
        "@x\nA\n-\nI\n", "@x\nAC\n+\nI\n", "@x\nA\n+\nI",
        "@x\nA\n+\n", "@x\nA\n", "@x\r\nA\n+\nI\n", std::string("@x\nA\n+\n\0\n", 10)
    };
    for (const auto &s : bad) rejects([&]() { fq::Frame x; x.append(s.data(), s.size()); x.finish(); });
    const std::string one = "@x\nA\n+\nI\n";
    fq::Frame exact(one.size()); exact.append(one.data(), one.size()); exact.finish();
    rejects([&]() { fq::Frame x(one.size() - 1); x.append(one.data(), one.size()); });
    fq::Frame failed(5);
    rejects([&]() { failed.append("@12345", 6); });
    rejects([&]() { failed.finish(); }); rejects([&]() { failed.append("\n", 1); });
    rejects([&]() { fq::Frame x(fq::RECORD_CAP, UINT64_MAX - 1); x.append(one.data(), one.size()); });
    rejects([&]() { fq::Frame x(fq::RECORD_CAP, 0, UINT64_MAX); x.append(one.data(), one.size()); x.take_prefix(1); });
    const std::string large = "@long\n" + std::string(400000, 'A') + "\n+\n" + std::string(400000, '@') + "\n";
    fq::Frame long_frame;
    for (size_t i = 0; i < large.size(); i += 37) long_frame.append(large.data() + i, std::min<size_t>(37, large.size() - i));
    long_frame.finish(); auto owned = long_frame.take_prefix(1);
    assert(owned.size == large.size() && !std::memcmp(owned.data, large.data(), large.size()));
    char *released = owned.release(); assert(!owned.data); std::free(released);
    auto name_check = [](const std::string &a, const std::string &b) {
        fq::Frame x, y; x.append(a.data(), a.size()); y.append(b.data(), b.size());
        x.finish(); y.finish(); auto p = x.take_prefix(1); auto q = y.take_prefix(1);
        return fq::same_name(reinterpret_cast<unsigned char *>(p.data), p.records[0],
                             reinterpret_cast<unsigned char *>(q.data), q.records[0]);
    };
    assert(name_check("@same a\nA\n+\n@\n", "@same\tb\nACG\n+\n+@@\n"));
    assert(!name_check("@same/1\nA\n+\n@\n", "@same/2\nA\n+\n@\n"));
    fq::LineOffsets overflow; overflow.scan(ptr(one), one.size());
    rejects([&]() { fq::split(overflow, UINT64_MAX - 1, true); });
    std::cout << "FASTQ_FRAME_TESTS_PASS exhaustive_cuts=" << (text.size() + 1) * (text.size() + 2) / 2 << '\n';
}
