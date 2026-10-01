#include "swbwa_sam_md5.h"
#include "swbwa_md5.hpp"
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <limits>

namespace {
struct Stream {
    Md5 digest;
    std::uint64_t bytes = 0, records = 0, headers = 0;
    unsigned char prefix[4] = {};
    unsigned prefix_size = 0;
    bool open = false, failed = false, leading = true, skipping = false;
    unsigned char last = 0;
} state;

int fail(int error) { state.failed = true; errno = error; return -1; }

int body(const unsigned char *data, size_t size)
{
    if (size > std::numeric_limits<std::uint64_t>::max() - state.bytes)
        return fail(EOVERFLOW);
    if (!size) return 0;
    state.digest.update(data, size);
    state.bytes += size;
    state.records += std::count(data, data + size, '\n');
    state.last = data[size - 1];
    return 0;
}

bool header(const unsigned char *p)
{
    return !std::memcmp(p, "@HD\t", 4) || !std::memcmp(p, "@SQ\t", 4) ||
           !std::memcmp(p, "@RG\t", 4) || !std::memcmp(p, "@PG\t", 4) ||
           !std::memcmp(p, "@CO\t", 4);
}
}

extern "C" int swbwa_sam_md5_open(void)
{
    if (state.open) return fail(EALREADY);
    state = Stream();
    state.open = true;
    return 0;
}

extern "C" int swbwa_sam_md5_update(const void *input, size_t size)
{
    const unsigned char *data = static_cast<const unsigned char *>(input);
    if (!state.open || state.failed || (!data && size)) return fail(EINVAL);
    if (!size) return 0;
    if (std::memchr(data, 0, size)) return fail(EILSEQ);
    while (size && state.leading) {
        if (state.skipping) {
            const unsigned char *newline = static_cast<const unsigned char *>(std::memchr(data, '\n', size));
            if (!newline) return 0;
            size_t used = newline - data + 1;
            data += used;
            size -= used;
            state.skipping = false;
            ++state.headers;
        } else {
            size_t used = std::min(size, size_t(4 - state.prefix_size));
            std::memcpy(state.prefix + state.prefix_size, data, used);
            state.prefix_size += used;
            data += used;
            size -= used;
            if (state.prefix_size != 4) return 0;
            if (header(state.prefix)) state.skipping = true;
            else {
                state.leading = false;
                if (body(state.prefix, 4)) return -1;
            }
            state.prefix_size = 0;
        }
    }
    return body(data, size);
}

extern "C" int swbwa_sam_md5_close(FILE *report)
{
    if (!state.open || state.failed || !report) return fail(EINVAL);
    if (state.skipping) return fail(EILSEQ);
    if (state.prefix_size && body(state.prefix, state.prefix_size)) return -1;
    if (!state.records || state.last != '\n') return fail(EILSEQ);
    const auto digest = state.digest.finish();
    char hex[33];
    const char *digits = "0123456789abcdef";
    for (unsigned i = 0; i < 16; ++i) {
        hex[2*i] = digits[digest[i] >> 4];
        hex[2*i+1] = digits[digest[i] & 15];
    }
    hex[32] = 0;
    state.open = false;
    int written = std::fprintf(report, "SAM_STREAM_MD5 {\"md5\": \"%s\", \"bytes\": %llu, "
        "\"records\": %llu, \"header_lines\": %llu, \"normalization\": "
        "\"leading header only; no sorting\", \"status\": \"COMPLETE\"}\n",
        hex, (unsigned long long)state.bytes, (unsigned long long)state.records,
        (unsigned long long)state.headers);
    return written < 0 || std::fflush(report) ? fail(EIO) : 0;
}

extern "C" void swbwa_sam_md5_abort(void) { state = Stream(); }

#ifdef SWBWA_SAM_MD5_STANDALONE
int main()
{
    static unsigned char buffer[1 << 20];
    if (swbwa_sam_md5_open()) return 1;
    size_t size;
    while ((size = std::fread(buffer, 1, sizeof(buffer), stdin)) != 0)
        if (swbwa_sam_md5_update(buffer, size)) { std::perror("SAM MD5"); return 1; }
    if (std::ferror(stdin) || swbwa_sam_md5_close(stdout)) {
        std::perror("SAM MD5"); return 1;
    }
    return 0;
}
#endif
