/* Native compute-node verifier: zlib validates every concatenated gzip member. */
#include "swbwa_md5.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

static void require_gzip_header(const std::string& path)
{
    if (path.size() < 3 || path.compare(path.size() - 3, 3, ".gz")) return;
    FILE *file = std::fopen(path.c_str(), "rb");
    if (!file) throw std::runtime_error("cannot open " + path);
    unsigned char magic[2] = {};
    bool valid = std::fread(magic, 1, 2, file) == 2 && magic[0] == 31 && magic[1] == 139;
    std::fclose(file);
    if (!valid) throw std::runtime_error("missing standard gzip header: " + path);
}

class Input {
    gzFile file;
    std::string path;
public:
    explicit Input(const std::string& name) : file(nullptr), path(name)
    {
        require_gzip_header(name);
        file = gzopen(name.c_str(), "rb");
        if (!file) throw std::runtime_error("cannot open " + name);
    }
    ~Input() { if (file) gzclose(file); }
    size_t read(char *buffer, unsigned size)
    {
        int n = gzread(file, buffer, size);
        if (n < 0) throw std::runtime_error("gzip read/checksum failed: " + path);
        return (size_t)n;
    }
    void finish()
    {
        int error = Z_OK;
        (void)gzerror(file, &error);
        int closed = gzclose(file);
        file = nullptr;
        if ((error != Z_OK && error != Z_STREAM_END) || closed != Z_OK)
            throw std::runtime_error("gzip checksum/close failed: " + path);
    }
};

static void report(Md5& hash, bool ordered, size_t files, size_t records,
                    size_t bytes, const std::string& first)
{
    auto digest = hash.finish();
    std::printf("FULL_SAM_COMPARE_PASS ordered=%d files=%zu records=%zu bytes=%zu md5=",
                ordered, files, records, bytes);
    for (auto byte : digest) std::printf("%02x", (unsigned)byte);
    std::printf(" first_file=%s\n", first.c_str());
}

static void check_ordered(const std::string& reference, const std::vector<std::string>& paths)
{
    Input expected(reference);
    char left[65536], right[65536];
    size_t records = 0, bytes = 0;
    Md5 hash;
    for (const auto& path : paths) {
        Input actual(path);
        size_t n;
        while ((n = actual.read(right, sizeof(right))) != 0) {
            size_t position = 0;
            while (position < n) {
                size_t part = expected.read(left + position, (unsigned)(n - position));
                if (!part) throw std::runtime_error("candidate has extra SAM bytes: " + path);
                position += part;
            }
            if (std::memcmp(left, right, n)) throw std::runtime_error("full SAM bytes differ: " + path);
            hash.update(right, n);
            bytes += n;
            records += std::count(right, right + n, '\n');
        }
        actual.finish();
    }
    if (expected.read(left, 1)) throw std::runtime_error("candidate is missing SAM bytes");
    expected.finish();
    report(hash, true, paths.size(), records, bytes, paths[0]);
}

static std::vector<std::string> read_lines(const std::vector<std::string>& paths)
{
    std::vector<std::string> lines;
    for (const auto& path : paths) {
        require_gzip_header(path);
        gzFile file = gzopen(path.c_str(), "rb");
        if (!file) throw std::runtime_error("cannot open " + path);
        char buffer[65536];
        std::string pending;
        int n;
        while ((n = gzread(file, buffer, sizeof(buffer))) > 0) {
            size_t begin = 0;
            for (int i = 0; i < n; ++i) if (buffer[i] == '\n') {
                pending.append(buffer + begin, (size_t)i + 1 - begin);
                lines.push_back(std::move(pending));
                pending.clear();
                begin = (size_t)i + 1;
            }
            pending.append(buffer + begin, (size_t)n - begin);
        }
        int error = 0;
        std::string message = gzerror(file, &error);
        int closed = gzclose(file);
        if (n < 0 || (error != Z_OK && error != Z_STREAM_END) || closed != Z_OK || !pending.empty())
            throw std::runtime_error("invalid/incomplete SAM or gzip " + path + ": " + message);
    }
    return lines;
}

int main(int argc, char **argv)
{
    try {
        bool unordered = argc > 1 && !std::strcmp(argv[1], "--unordered");
        bool each = argc > 1 && !std::strcmp(argv[1], "--each");
        int first = unordered || each ? 2 : 1;
        if (argc < first + 2) throw std::runtime_error("usage: check_gzip_sam [--unordered|--each] reference candidate...");
        std::vector<std::string> paths;
        for (int i = first + 1; i < argc; ++i) paths.emplace_back(argv[i]);
        if (!unordered) {
            if (each) for (const auto& path : paths) check_ordered(argv[first], {path});
            else check_ordered(argv[first], paths);
            return 0;
        }
        auto expected = read_lines({argv[first]});
        if (unordered) std::sort(expected.begin(), expected.end());
        for (size_t index = 0; index < (each ? paths.size() : 1); ++index) {
            const std::vector<std::string> current = each ? std::vector<std::string>{paths[index]} : paths;
            auto actual = read_lines(current);
            if (unordered) std::sort(actual.begin(), actual.end());
            if (actual != expected) throw std::runtime_error("full SAM bytes/records differ: " + current[0]);
            Md5 hash;
            size_t bytes = 0;
            for (const auto& line : actual) {
                hash.update(line.data(), line.size());
                bytes += line.size();
            }
            report(hash, false, current.size(), actual.size(), bytes, current[0]);
        }
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "FULL_SAM_COMPARE_FAIL %s\n", error.what());
        return 1;
    }
}
