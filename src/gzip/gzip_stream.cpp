#include "gzip_stream.h"
#include "rapid_reader.h"
#include "byte_source.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

namespace {
double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
void validate(size_t chunk, unsigned workers) {
    if (chunk < (64U << 10) || chunk > (64U << 20) || !workers || workers > 384)
        throw std::runtime_error("invalid gzip chunk size or worker count");
}
void launch(std::vector<swbwa_codec_job> &jobs, unsigned workers, SwbwaGzipMetrics &m) {
    if (jobs.empty()) return;
    swbwa_codec_batch batch{jobs.data(), (unsigned)jobs.size(), workers, 0, 0};
    double t = now();
    swbwa_codec_launch(&batch);
    m.launch_s += now() - t;
    m.kernel_s += batch.run_s; m.queue_s += batch.wait_s;
    m.jobs += jobs.size();
    for (const auto &job : jobs)
        if (job.status != Z_OK)
            throw std::runtime_error("parallel gzip job failed, status=" + std::to_string(job.status));
}
}

struct SwbwaGzipReader::Impl {
    std::unique_ptr<SwbwaRapidReader> rapid;
    std::unique_ptr<SwbwaByteSource> serial;
    std::vector<unsigned char> input;
    z_stream stream{};
    bool initialized = false;
    bool ended = false;
    bool member_end = false;
    bool views = false;
    ~Impl() {
        if (initialized) inflateEnd(&stream);
    }
    void refill() {
        if (stream.avail_in) return;
        stream.avail_in = (uInt)serial->read(input.data(), input.size());
        stream.next_in = input.data();
    }
    size_t read_serial(void *out, size_t bytes) {
        if (ended) return 0;
        const uInt capacity = (uInt)std::min<size_t>(bytes, INT_MAX);
        stream.next_out = static_cast<Bytef *>(out);
        stream.avail_out = capacity;
        while (stream.avail_out && !ended) {
            refill();
            if (member_end) {
                if (!stream.avail_in) { ended = true; break; }
                if (*stream.next_in == 0) {
                    do {
                        while (stream.avail_in) {
                            if (*stream.next_in++ != 0)
                                throw std::runtime_error("nonzero bytes after gzip padding");
                            --stream.avail_in;
                        }
                        refill();
                    } while (stream.avail_in);
                    ended = true;
                    break;
                }
                Bytef *next_in = stream.next_in, *next_out = stream.next_out;
                uInt available_in = stream.avail_in, available_out = stream.avail_out;
                if (inflateReset(&stream) != Z_OK) throw std::runtime_error("gzip reset failed");
                stream.next_in = next_in; stream.avail_in = available_in;
                stream.next_out = next_out; stream.avail_out = available_out;
                member_end = false;
            }
            uInt before_in = stream.avail_in, before_out = stream.avail_out;
            int result = inflate(&stream, Z_NO_FLUSH);
            if (result == Z_STREAM_END) {
                member_end = true;
            } else if (result != Z_OK ||
                       (before_in == stream.avail_in && before_out == stream.avail_out)) {
                throw std::runtime_error("invalid or truncated gzip stream");
            }
        }
        return capacity - stream.avail_out;
    }
};

SwbwaGzipReader::SwbwaGzipReader(const char *path, size_t chunk, unsigned workers,
                             const std::string &backend,
                             const std::vector<unsigned char> *memory, bool views) : impl_(new Impl) {
    validate(chunk, workers);
    if (backend != "auto" && backend != "rapid" && backend != "zlib")
        throw std::runtime_error("gzip reader must be auto, rapid or zlib");
    bool rapid = backend == "rapid";
#ifdef SWBWA_CODEC_RAPID
    rapid = rapid || backend == "auto";
#else
    if (rapid) throw std::runtime_error("build with GZIP_DECODER=rapid for parallel input");
#endif
    if (rapid) {
        impl_->views = views;
        impl_->rapid.reset(new SwbwaRapidReader(path, chunk, workers, metrics_, memory, views));
    }
    else {
        impl_->serial.reset(new SwbwaByteSource(path, memory));
        metrics_.input_bytes = impl_->serial->size();
        metrics_.memory_source = memory != nullptr;
        impl_->input.resize(1U << 20);
        if (inflateInit2(&impl_->stream, 31) != Z_OK) throw std::runtime_error("gzip initialize failed");
        impl_->initialized = true;
    }
}
SwbwaGzipReader::~SwbwaGzipReader() = default;
bool SwbwaGzipReader::supports_views() const { return impl_->views; }
bool SwbwaGzipReader::peek(size_t maximum, SwbwaDecodedView &view) {
    if (!supports_views()) throw std::logic_error("gzip reader does not support views");
    const double start = now();
    bool present = impl_->rapid->peek(maximum, view);
    metrics_.read_s += now() - start;
    return present;
}
void SwbwaGzipReader::consume(size_t bytes) {
    if (!supports_views()) throw std::logic_error("gzip reader does not support views");
    impl_->rapid->consume(bytes);
    metrics_.output_bytes += bytes;
}
bool SwbwaGzipReader::eof() const {
    return impl_->rapid ? impl_->rapid->eof() : impl_->ended;
}
size_t SwbwaGzipReader::read(void *destination, size_t bytes, SwbwaLineIndex *lines) {
    if (!bytes) return 0;
    const double start = now();
    size_t copied;
    if (impl_->rapid) copied = impl_->rapid->read(destination, bytes, lines);
    else {
        copied = impl_->read_serial(destination, bytes);
        metrics_.serial_s += now() - start;
    }
    if (lines && !impl_->rapid) {
        const double index_start = now();
        lines->append(static_cast<unsigned char *>(destination), copied);
        metrics_.line_index_s += now() - index_start;
    }
    metrics_.read_s += now() - start;
    metrics_.output_bytes += copied;
    return copied;
}

SwbwaGzipWriter::SwbwaGzipWriter(const char *path, size_t chunk, unsigned workers, int level)
    : path_(path), chunk_bytes_(chunk), workers_(workers), level_(level) {
    validate(chunk, workers);
    if (level < 1 || level > 9) throw std::runtime_error("invalid gzip compression level");
    batch_chunks_ = std::max<size_t>(1, std::min<size_t>(workers, (64U << 20) / chunk));
    int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd < 0) throw std::runtime_error("gzip output already exists or cannot be created: " + path_);
    file_ = fdopen(fd, "wb");
    if (!file_) { close(fd); unlink(path); throw std::runtime_error("gzip fdopen failed"); }
    const unsigned char header[10] = {31,139,8,0,0,0,0,0,0,3};
    try { write(header, sizeof(header)); }
    catch (...) { fclose(file_); file_ = nullptr; unlink(path); throw; }
}
SwbwaGzipWriter::~SwbwaGzipWriter() {
    if (file_) fclose(file_);
    if (!complete_) unlink(path_.c_str());
}
void SwbwaGzipWriter::write(const void *data, size_t bytes) {
    double t = now();
    if (bytes && fwrite(data, 1, bytes, file_) != bytes) throw std::runtime_error("gzip output write failed");
    metrics_.io_s += now() - t;
    metrics_.output_bytes += bytes;
}
void SwbwaGzipWriter::append(const void *data, size_t bytes) {
    if (!file_) throw std::runtime_error("gzip append after close");
    const auto *p = static_cast<const unsigned char *>(data);
    while (bytes) {
        if (input_.empty() || input_.back().size() == chunk_bytes_) input_.emplace_back();
        auto &b = input_.back();
        size_t n = std::min(bytes, chunk_bytes_ - b.size());
        b.insert(b.end(), p, p + n);
        p += n; bytes -= n; metrics_.input_bytes += n;
        if (input_.size() == batch_chunks_ && b.size() == chunk_bytes_) flush();
    }
}
void SwbwaGzipWriter::flush() {
    std::vector<std::vector<unsigned char>> output(input_.size()), dictionaries(input_.size());
    std::vector<swbwa_codec_job> jobs;
    for (size_t i = 0; i < input_.size(); ++i) {
        const auto &b = input_[i];
        dictionaries[i] = previous_;
        // Append (not replace) history when the final block is shorter than 32 KiB.
        if (b.size() >= 32768) previous_.assign(b.end() - 32768, b.end());
        else {
            previous_.insert(previous_.end(), b.begin(), b.end());
            if (previous_.size() > 32768) previous_.erase(previous_.begin(), previous_.end() - 32768);
        }
        output[i].resize(compressBound(b.size()) + 64);
        swbwa_codec_job job{};
        job.input = b.data(); job.input_size = b.size();
        job.output = output[i].data(); job.output_capacity = output[i].size();
        job.dictionary = dictionaries[i].data(); job.dictionary_size = dictionaries[i].size();
        job.level = level_; job.operation = SWBWA_CODEC_DEFLATE;
        jobs.push_back(job);
    }
    launch(jobs, workers_, metrics_);
    for (const auto &job : jobs) {
        write(job.output, job.output_size);
        crc_ = (uint32_t)crc32_combine(crc_, job.crc, (z_off_t)job.input_size);
    }
    input_.clear();
}
void SwbwaGzipWriter::finish() {
    if (complete_) return;
    flush();
    unsigned char trailer[10] = {3,0};
    for (unsigned i = 0; i < 4; ++i) {
        trailer[2+i] = (unsigned char)(crc_ >> (8*i));
        trailer[6+i] = (unsigned char)((uint32_t)metrics_.input_bytes >> (8*i));
    }
    write(trailer, sizeof(trailer));
    FILE *p = file_; file_ = nullptr;
    if (fclose(p)) throw std::runtime_error("gzip close failed");
    complete_ = true;
}
