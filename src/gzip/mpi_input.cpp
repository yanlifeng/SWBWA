#include "swbwa_mpi_gzip.h"
#include "swbwa_config.h"
#include "mpi_decoder.h"
#include "fastq_frame.h"
#include "line_index.h"
#include "swbwa_host_workers.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace fq = swbwa_fastq;
typedef unsigned long long U64;
char last_error[512] = "";

void check(int code)
{
    if (code != MPI_SUCCESS) {
        char message[MPI_MAX_ERROR_STRING]; int n = 0;
        MPI_Error_string(code, message, &n);
        throw std::runtime_error(std::string("MPI gzip: ") + std::string(message, size_t(n)));
    }
}
void thread_call(int code)
{
    if (code) throw std::runtime_error(std::string("MPI gzip pthread: ") + std::strerror(code));
}
uint64_t add(uint64_t a, uint64_t b)
{
    if (b > UINT64_MAX - a) throw std::overflow_error("MPI gzip ordinal/byte overflow");
    return a + b;
}
uint64_t multiply(uint64_t a, uint64_t b)
{
    if (b && a > UINT64_MAX / b) throw std::overflow_error("MPI gzip ordinal/byte overflow");
    return a * b;
}
void thread_check()
{
    int level = 0, main_thread = 0;
    check(MPI_Query_thread(&level)); check(MPI_Is_thread_main(&main_thread));
    if (!main_thread && level < MPI_THREAD_MULTIPLE)
        throw std::runtime_error("MPI gzip worker thread requires MPI_THREAD_MULTIPLE");
}
void put(std::vector<unsigned char> &out, uint64_t v)
{
    for (unsigned i = 0; i < 8; ++i) { out.push_back(static_cast<unsigned char>(v)); v >>= 8; }
}
uint64_t get(const unsigned char *&p, const unsigned char *end)
{
    if (end - p < 8) throw std::runtime_error("truncated MPI FASTQ packet");
    uint64_t n = 0;
    for (unsigned i = 0; i < 8; ++i) n |= uint64_t(*p++) << (8 * i);
    return n;
}
struct Piece {
    fq::OwnedChunk buffer;
    SwbwaLineIndex lines;
    uint64_t records = 0;
};
struct Chunk {
    std::map<uint64_t, Piece> pieces;
    uint64_t records = 0, bytes = 0, charge = 0;
};
struct Stream {
    swbwa_mpi_decoder *decoder = NULL;
    uint64_t bytes = 0, lfs = 0, records = 0, carry_size = 0, queued = 0;
    uint64_t epochs = 0, received_records = 0, peak_queued = 0, boundary_records = 0;
    int carry_owner = -1;
    bool eof = false, after_lf = true;
    std::vector<unsigned char> carry;
    std::map<uint64_t, Chunk> chunks;
};
struct Local {
    uint64_t first = 0;
    std::vector<unsigned char> crossing;
    SwbwaLineIndex crossing_lines, body_lines;
    size_t body_begin = 0;
};

struct ReadyChunk {
    swbwa_mpi_gzip_batch batch = {};
    ~ReadyChunk()
    {
        std::free(batch.data[0]); std::free(batch.data[1]); std::free(batch.partitions);
    }
    void *buffer(unsigned i) const
    {
        return i < 2 ? static_cast<void *>(batch.data[i]) : static_cast<void *>(batch.partitions);
    }
};

struct ReadyDescriptor {
    /* ID, ordinal, records, mate lengths, partition count, three displacements. */
    std::array<U64, 9> values = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
};

unsigned ready_chunk_limit()
{
    const char *value = std::getenv("SWBWA_GZIP_READY_CHUNKS");
    if (!value) return 4;
    char *end = NULL;
    errno = 0;
    const unsigned long n = std::strtoul(value, &end, 10);
    if (errno || end == value || *end || n < 1 || n > 16)
        throw std::invalid_argument("SWBWA_GZIP_READY_CHUNKS must be 1..16");
    return unsigned(n);
}

fq::Split split_index(const unsigned char *data, const SwbwaLineIndex &index,
                      uint64_t before, bool after_lf)
{
    const bool boundary = !(before % 4) && after_lf;
    fq::Split s = {0, 0, 0, 0, false, before / 4};
    const size_t first = 4 - size_t(before % 4);
    if (first > index.lines()) {
        if (!boundary) s.prefix_end = s.body_begin = s.body_end = s.suffix_begin = index.bytes();
        return s;
    }
    if (!boundary) {
        s.prefix_end = s.body_begin = index.end_of_line(data, first);
        s.prefix_completes = true;
    }
    const size_t last = first + ((index.lines() - first) / 4) * 4;
    s.body_end = s.suffix_begin = index.end_of_line(data, last);
    return s;
}

struct Copy { char *to; const unsigned char *from; size_t size; };
void copy_slices(void *opaque, int worker, int workers)
{
    Copy &copy = *static_cast<Copy *>(opaque);
    const size_t begin = copy.size / size_t(workers) * size_t(worker);
    const size_t end = worker + 1 == workers ? copy.size :
        copy.size / size_t(workers) * size_t(worker + 1);
    if (end > begin) std::memcpy(copy.to + begin, copy.from + begin, end - begin);
}
void copy_bytes(char *to, const unsigned char *from, size_t size)
{
    Copy copy = {to, from, size};
    if (size >= (1U << 20)) swbwa_host_workers_run_ready(copy_slices, &copy);
    else if (size) std::memcpy(to, from, size);
}
}

struct swbwa_mpi_gzip_reader {
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0, ranks = 0, mates = 1, debug = 0;
    unsigned workers = 0;
    size_t target = 0;
    uint64_t k = 0, wave = 0, queue_limit = 0, emitted = 0, emitted_chunks = 0, empty = 0;
    uint64_t emitted_bytes[2] = {0, 0};
    bool done = false, decoded_done = false;
    bool dynamic = SWBWA_MPI_INPUT_MODE == SWBWA_MPI_INPUT_DYNAMIC;
    unsigned ready_limit = 0;
    MPI_Win ticket_window = MPI_WIN_NULL, data_window = MPI_WIN_NULL;
    U64 ticket = 0;
    uint64_t ready_begin = 0, ready_count = 0, ready_windows = 0;
    uint64_t ready_peak_bytes = 0, local_claims = 0, remote_claims = 0, fetched_bytes = 0;
    double refill_s = 0, refill_wait_s = 0, drain_wait_s = 0, claim_s = 0, fetch_s = 0;
    double pack_s = 0, local_copy_s = 0;
    pthread_mutex_t pending_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t completed = PTHREAD_COND_INITIALIZER;
    std::deque<int64_t> pending;
    fq::OwnedChunk ready_payload;
    std::vector<ReadyDescriptor> descriptors;
    double decode_s = 0, frame_s = 0, exchange_s = 0, emit_s = 0, next_s = 0;
    double frame_mpi_s = 0, exchange_mpi_s = 0, ready_wait_s = 0, eof_wait_s = 0;
    uint64_t local_payload_bytes = 0, remote_payload_bytes = 0, assembled_copy_bytes = 0;
    uint64_t mate_same_offset = 0, mate_corrected_offset = 0;
    char error[512] = "";
    Stream streams[2];

    void advance(int mate);
    void exchange(Stream &, const Local &, const unsigned char *);
    void pack(Stream &, std::vector<std::vector<unsigned char> > &, uint64_t,
              const unsigned char *, const SwbwaLineIndex &);
    void insert(Stream &, uint64_t, Piece &&);
    void receive(Stream &, const unsigned char *, size_t);
    fq::OwnedChunk assemble(Stream &, uint64_t id, uint64_t first, uint64_t count, SwbwaLineIndex &);
    void partitions(const fq::OwnedChunk *, const SwbwaLineIndex *, uint64_t,
                    swbwa_fastq_partition_t *);
    int next(swbwa_mpi_gzip_batch &);
    int next_static(swbwa_mpi_gzip_batch &);
    int next_dynamic(swbwa_mpi_gzip_batch &);
    bool refill();
    void release_ready();
    ~swbwa_mpi_gzip_reader()
    {
        pthread_cond_destroy(&completed);
        pthread_mutex_destroy(&pending_mutex);
    }
};

namespace {
int abort_reader(swbwa_mpi_gzip_reader *r, const char *message)
{
    std::snprintf(last_error, sizeof(last_error), "%s", message);
    if (r) std::snprintf(r->error, sizeof(r->error), "%s", message);
    std::fprintf(stderr, "[MPI gzip rank=%d] fatal: %s\n", r ? r->rank : -1, message);
    std::fflush(stderr);
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Finalized(&finalized);
    if (initialized && !finalized) {
        MPI_Abort(r && r->comm != MPI_COMM_NULL ? r->comm : MPI_COMM_WORLD, 1);
        std::abort(); // MPI_Abort must not return into a partially failed collective.
    }
    return -1;
}
}

void swbwa_mpi_gzip_reader::advance(int mate)
{
    Stream &s = streams[mate];
    if (s.eof) throw std::logic_error("advance after decoder EOF");
    const double t0 = MPI_Wtime();
    swbwa_mpi_decoded_batch batch = {};
    const int status = swbwa_mpi_decoder_next(s.decoder, &batch);
    if (status < 0) throw std::runtime_error(swbwa_mpi_decoder_error(s.decoder));
    if (!status) throw std::runtime_error("decoder EOF without committed final epoch");
    decode_s += MPI_Wtime() - t0;
    const double t1 = MPI_Wtime();
    SwbwaLineIndex lines;
    for (size_t i = 0; i < batch.count; ++i) {
        const swbwa_mpi_decoded_extent &extent = batch.extents[i];
        lines.append(extent.data, extent.size, extent.line_counts);
    }
    if (lines.bytes() != batch.size) throw std::runtime_error("decoded newline metadata size mismatch");
    const std::array<U64, 7> own = {{batch.decoded_offset, batch.size, lines.lines(),
        batch.size ? U64(batch.bytes[batch.size - 1] == '\n') : 0,
        batch.global_begin, batch.global_decoded_end, U64(batch.eof != 0)}};
    std::vector<std::array<U64, 7> > metadata(size_t(ranks), std::array<U64, 7>{{}});
    const double gather_begin = MPI_Wtime();
    check(MPI_Allgather(own.data(), 7, MPI_UNSIGNED_LONG_LONG, metadata.data(), 7,
                        MPI_UNSIGNED_LONG_LONG, comm));
    frame_mpi_s += MPI_Wtime() - gather_begin;
    std::vector<int> order;
    for (int i = 0; i < ranks; ++i) {
        if (metadata[i][4] != s.bytes || metadata[i][5] != batch.global_decoded_end ||
            metadata[i][6] != U64(batch.eof != 0))
            throw std::runtime_error("inconsistent decoder epoch metadata");
        if (metadata[i][1]) order.push_back(i);
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) { return metadata[a][0] < metadata[b][0]; });
    uint64_t offset = s.bytes, lf = s.lfs, before = 0;
    bool after = s.after_lf, local_after = false;
    for (int owner : order) {
        if (metadata[owner][0] != offset) throw std::runtime_error("noncontiguous decoded extents");
        if (owner == rank) { before = lf; local_after = after; }
        offset = add(offset, metadata[owner][1]); lf = add(lf, metadata[owner][2]);
        after = metadata[owner][3] != 0;
    }
    if (offset != batch.global_decoded_end) throw std::runtime_error("decoded extent end mismatch");
    fq::Split split = {};
    Local local;
    if (batch.size) {
        split = split_index(batch.bytes, lines, before, local_after);
        local.first = split.first_record;
        local.body_begin = split.body_begin;
        local.body_lines.append_slice(lines, batch.bytes + split.body_begin,
                                      split.body_begin, split.body_end - split.body_begin);
    }

    // Only a bounded, unfinished record follows file order. Interior records
    // stay as bulk byte ranges plus newline counts produced by the CPE decoder.
    int previous = s.carry_owner;
    uint64_t carry_size = s.carry_size;
    for (int owner : order) {
        if (carry_size && previous != owner) {
            if (rank == previous) {
                if (s.carry.size() != carry_size) throw std::logic_error("missing boundary carry");
                const double start = MPI_Wtime();
                check(MPI_Send(s.carry.data(), int(carry_size), MPI_BYTE, owner, 41, comm));
                frame_mpi_s += MPI_Wtime() - start;
                s.carry.clear();
            }
            if (rank == owner) {
                s.carry.resize(size_t(carry_size));
                const double start = MPI_Wtime();
                check(MPI_Recv(s.carry.data(), int(carry_size), MPI_BYTE, previous, 41, comm, MPI_STATUS_IGNORE));
                frame_mpi_s += MPI_Wtime() - start;
            }
        }
        U64 outgoing = 0;
        if (rank == owner) {
            if (!carry_size) s.carry.clear();
            if (split.prefix_end) {
                if (split.prefix_end > fq::RECORD_CAP - s.carry.size())
                    throw std::runtime_error("FASTQ boundary record exceeds cap");
                s.carry.insert(s.carry.end(), batch.bytes, batch.bytes + split.prefix_end);
                if (split.prefix_completes) {
                    local.crossing.swap(s.carry);
                    fq::LineOffsets crossing_lines;
                    crossing_lines.scan(local.crossing.data(), local.crossing.size());
                    if (fq::parse(local.crossing.data(), crossing_lines, 0, local.crossing.size()).size() != 1)
                        throw std::runtime_error("FASTQ boundary carry does not contain one record");
                    local.crossing_lines.append(local.crossing.data(), local.crossing.size());
                }
            } else if (!s.carry.empty()) {
                throw std::runtime_error("FASTQ phase/carry mismatch");
            }
            if (!split.prefix_end || split.prefix_completes) {
                const size_t n = batch.size - split.suffix_begin;
                if (n > fq::RECORD_CAP) throw std::runtime_error("FASTQ suffix exceeds cap");
                s.carry.assign(batch.bytes + split.suffix_begin, batch.bytes + batch.size);
            }
            outgoing = s.carry.size();
        }
        const double broadcast_begin = MPI_Wtime();
        check(MPI_Bcast(&outgoing, 1, MPI_UNSIGNED_LONG_LONG, owner, comm));
        frame_mpi_s += MPI_Wtime() - broadcast_begin;
        if (outgoing > fq::RECORD_CAP) throw std::runtime_error("FASTQ carry exceeds cap");
        carry_size = outgoing; previous = owner;
    }
    s.bytes = offset; s.lfs = lf; s.records = lf / 4; s.after_lf = after;
    s.carry_owner = previous; s.carry_size = carry_size; s.eof = batch.eof != 0; ++s.epochs;
    s.boundary_records += !local.crossing.empty();
    if (s.eof && (s.carry_size || s.lfs % 4))
        throw std::runtime_error("truncated FASTQ: four LF-terminated lines required at EOF");
    if (!k && mate == 0 && (s.records || s.eof)) {
        const uint64_t complete_bytes = s.bytes - s.carry_size;
        const uint64_t mean = s.records ? complete_bytes / s.records + (complete_bytes % s.records != 0) : 1;
        k = std::max<uint64_t>(1, target / std::max<uint64_t>(1, mean));
        k = std::min<uint64_t>(k, SWBWA_MAX_READS_PER_BATCH / mates);
        if (debug) std::fprintf(stderr, "[MPI gzip rank=%d] records_per_chunk=%llu target_R1_bytes=%zu\n",
                                rank, U64(k), target);
    }
    frame_s += MPI_Wtime() - t1;
    if (k) exchange(s, local, batch.bytes);
    else if (local.body_lines.lines() || !local.crossing.empty())
        throw std::logic_error("records exist before chunk estimate");
}

void swbwa_mpi_gzip_reader::insert(Stream &s, uint64_t first, Piece &&piece)
{
    const uint64_t count = piece.records, bytes = piece.buffer.size, id = first / k;
    if (!count || count > k - first % k || id % uint64_t(ranks) != uint64_t(rank) ||
        id / uint64_t(ranks) < wave || piece.lines.bytes() != bytes ||
        piece.lines.lines() != multiply(count, 4))
        throw std::runtime_error("invalid MPI FASTQ ordinal piece");
    const uint64_t charge = add(256, add(bytes, multiply(piece.lines.spans(), 2 * sizeof(size_t))));
    if (charge > queue_limit - s.queued) throw std::runtime_error("MPI FASTQ queue budget exceeded");
    Chunk &chunk = s.chunks[id];
    auto next = chunk.pieces.lower_bound(first);
    if (next != chunk.pieces.end() && add(first, count) > next->first)
        throw std::runtime_error("overlapping MPI FASTQ records");
    if (next != chunk.pieces.begin()) {
        auto previous = next; --previous;
        if (add(previous->first, previous->second.records) > first)
            throw std::runtime_error("overlapping MPI FASTQ records");
    }
    chunk.pieces.emplace_hint(next, first, std::move(piece));
    chunk.records = add(chunk.records, count); chunk.bytes = add(chunk.bytes, bytes);
    chunk.charge += charge; s.queued += charge; s.received_records = add(s.received_records, count);
    s.peak_queued = std::max(s.peak_queued, s.queued);
}

void swbwa_mpi_gzip_reader::pack(Stream &s, std::vector<std::vector<unsigned char> > &send,
    uint64_t first, const unsigned char *bytes, const SwbwaLineIndex &lines)
{
    if (lines.lines() % 4) throw std::runtime_error("incomplete FASTQ bulk interval");
    size_t begin = 0, byte_begin = 0;
    while (begin < lines.lines() / 4) {
        const uint64_t ordinal = add(first, begin), id = ordinal / k;
        const size_t count = size_t(std::min<uint64_t>(lines.lines() / 4 - begin, k - ordinal % k));
        const size_t byte_end = lines.end_of_line(bytes, (begin + count) * 4);
        SwbwaLineIndex slice;
        slice.append_slice(lines, bytes + byte_begin, byte_begin, byte_end - byte_begin);
        const int owner = int(id % uint64_t(ranks));
        if (owner == rank) {
            Piece piece; piece.records = count; piece.lines = std::move(slice);
            piece.buffer.size = byte_end - byte_begin;
            piece.buffer.data = static_cast<char *>(std::malloc(piece.buffer.size + 1));
            if (!piece.buffer.data) throw std::bad_alloc();
            copy_bytes(piece.buffer.data, bytes + byte_begin, piece.buffer.size);
            piece.buffer.data[piece.buffer.size] = 0;
            local_payload_bytes += piece.buffer.size;
            insert(s, ordinal, std::move(piece));
            begin += count; byte_begin = byte_end;
            continue;
        }
        std::vector<unsigned char> &out = send[size_t(owner)];
        const uint64_t packet_size = add(32, add(multiply(slice.spans(), 16), byte_end - byte_begin));
        if (add(out.size(), packet_size) > std::min<uint64_t>(INT_MAX, queue_limit))
            throw std::runtime_error("MPI FASTQ send budget exceeded");
        put(out, ordinal); put(out, count); put(out, byte_end - byte_begin);
        put(out, slice.spans());
        for (size_t i = 0; i < slice.spans(); ++i) {
            put(out, slice.span_bytes(i)); put(out, slice.span_lines(i));
        }
        out.insert(out.end(), bytes + byte_begin, bytes + byte_end);
        remote_payload_bytes += byte_end - byte_begin;
        begin += count; byte_begin = byte_end;
    }
}

void swbwa_mpi_gzip_reader::receive(Stream &s, const unsigned char *data, size_t size)
{
    const unsigned char *p = data, *end = data + size;
    while (p < end) {
        const uint64_t first = get(p, end), count = get(p, end), bytes = get(p, end), id = first / k;
        const uint64_t spans = get(p, end);
        if (!count || count > k - first % k || id % uint64_t(ranks) != uint64_t(rank) ||
            id / uint64_t(ranks) < wave || !spans || spans > uint64_t(end - p) / 16 || bytes >= SIZE_MAX)
            throw std::runtime_error("invalid MPI FASTQ ordinal packet");
        const uint64_t charge = add(256, add(bytes, multiply(spans, 2 * sizeof(size_t))));
        if (charge > queue_limit - s.queued) throw std::runtime_error("MPI FASTQ queue budget exceeded");
        Piece piece;
        piece.records = count;
        for (uint64_t i = 0; i < spans; ++i) {
            const uint64_t n = get(p, end), lfs = get(p, end);
            piece.lines.append_span(size_t(n), size_t(lfs));
        }
        if (piece.lines.bytes() != bytes || bytes > uint64_t(end - p))
            throw std::runtime_error("invalid MPI FASTQ payload size");
        piece.buffer.size = size_t(bytes);
        piece.buffer.data = static_cast<char *>(std::malloc(size_t(bytes) + 1));
        if (!piece.buffer.data) throw std::bad_alloc();
        copy_bytes(piece.buffer.data, p, size_t(bytes));
        piece.buffer.data[bytes] = 0; p += size_t(bytes);
        insert(s, first, std::move(piece));
    }
}

void swbwa_mpi_gzip_reader::exchange(Stream &s, const Local &local, const unsigned char *body)
{
    const double begin = MPI_Wtime();
    std::vector<std::vector<unsigned char> > send(static_cast<size_t>(ranks));
    pack(s, send, local.first, local.crossing.data(), local.crossing_lines);
    pack(s, send, add(local.first, !local.crossing.empty()),
         body ? body + local.body_begin : NULL, local.body_lines);
    if (ranks == 1) { exchange_s += MPI_Wtime() - begin; return; }
    std::vector<int> sc(size_t(ranks), 0), rc(size_t(ranks), 0), sd(size_t(ranks), 0), rd(size_t(ranks), 0);
    uint64_t send_size = 0, recv_size = 0;
    for (int i = 0; i < ranks; ++i) {
        sd[i] = int(send_size); sc[i] = int(send[i].size()); send_size = add(send_size, send[i].size());
        if (send_size > std::min<uint64_t>(INT_MAX, queue_limit))
            throw std::runtime_error("MPI FASTQ total send budget exceeded");
    }
    double mpi_begin = MPI_Wtime();
    check(MPI_Alltoall(sc.data(), 1, MPI_INT, rc.data(), 1, MPI_INT, comm));
    exchange_mpi_s += MPI_Wtime() - mpi_begin;
    for (int i = 0; i < ranks; ++i) {
        if (rc[i] < 0) throw std::runtime_error("negative MPI FASTQ receive size");
        rd[i] = int(recv_size); recv_size = add(recv_size, uint64_t(rc[i]));
        if (recv_size > std::min<uint64_t>(INT_MAX, queue_limit - s.queued))
            throw std::runtime_error("MPI FASTQ receive/queue budget exceeded");
    }
    std::vector<unsigned char> packed(size_t(send_size), 0), incoming(size_t(recv_size), 0);
    for (int i = 0; i < ranks; ++i) {
        if (sc[i]) std::memcpy(packed.data() + sd[i], send[i].data(), size_t(sc[i]));
        std::vector<unsigned char>().swap(send[i]);
    }
    unsigned char dummy = 0;
    mpi_begin = MPI_Wtime();
    check(MPI_Alltoallv(packed.empty() ? &dummy : packed.data(), sc.data(), sd.data(), MPI_BYTE,
                        incoming.empty() ? &dummy : incoming.data(), rc.data(), rd.data(), MPI_BYTE, comm));
    exchange_mpi_s += MPI_Wtime() - mpi_begin;
    for (int i = 0; i < ranks; ++i)
        if (rc[i]) receive(s, incoming.data() + rd[i], size_t(rc[i]));
    exchange_s += MPI_Wtime() - begin;
}

fq::OwnedChunk swbwa_mpi_gzip_reader::assemble(Stream &s, uint64_t id, uint64_t first,
                                              uint64_t count, SwbwaLineIndex &index)
{
    fq::OwnedChunk out; out.first_record = first;
    auto entry = s.chunks.find(id);
    if (!count) {
        if (entry != s.chunks.end()) throw std::runtime_error("unexpected records in empty MPI chunk");
        return out;
    }
    if (entry == s.chunks.end() || entry->second.records != count)
        throw std::runtime_error("missing MPI FASTQ records in complete chunk");
    Chunk &chunk = entry->second;
    if (chunk.bytes > uint64_t(SWBWA_CPE_FORMAT_BUFFER_BYTES) || chunk.bytes >= SIZE_MAX)
        throw std::runtime_error("MPI FASTQ chunk exceeds CPE formatting buffer");
    if (chunk.pieces.size() == 1 && chunk.pieces.begin()->first == first) {
        Piece &piece = chunk.pieces.begin()->second;
        out = std::move(piece.buffer); out.first_record = first;
        index = std::move(piece.lines);
        s.queued -= chunk.charge; s.chunks.erase(entry);
        return out;
    }
    out.size = size_t(chunk.bytes);
    out.data = static_cast<char *>(std::malloc(out.size + 1));
    if (!out.data) throw std::bad_alloc();
    size_t offset = 0; uint64_t ordinal = first;
    for (const auto &part : chunk.pieces) {
        if (part.first != ordinal) throw std::runtime_error("gap in MPI FASTQ record ordinals");
        const Piece &piece = part.second;
        copy_bytes(out.data + offset, reinterpret_cast<const unsigned char *>(piece.buffer.data), piece.buffer.size);
        assembled_copy_bytes += piece.buffer.size;
        index.append_index(piece.lines);
        offset += piece.buffer.size; ordinal = add(ordinal, piece.records);
    }
    if (offset != out.size || ordinal != add(first, count)) throw std::logic_error("MPI FASTQ assembly mismatch");
    out.data[out.size] = 0;
    s.queued -= chunk.charge; s.chunks.erase(entry);
    return out;
}

void swbwa_mpi_gzip_reader::partitions(const fq::OwnedChunk *buffers,
    const SwbwaLineIndex *lines, uint64_t count, swbwa_fastq_partition_t *out)
{
    uint64_t previous = 0;
    size_t offsets[2] = {0, 0};
    for (unsigned i = 0; i < workers; ++i) {
        uint64_t record = count;
        size_t ends[2] = {buffers[0].size, buffers[1].size};
        if (i + 1 != workers) {
            const size_t nominal = buffers[0].size / workers * (i + 1);
            record = (lines[0].lines_before(
                reinterpret_cast<const unsigned char *>(buffers[0].data), nominal) + 3) / 4;
            record = std::min(count, record);
            ends[0] = record ? lines[0].end_of_line(
                reinterpret_cast<const unsigned char *>(buffers[0].data), size_t(record * 4)) : 0;
            if (mates == 2) {
                const unsigned char *r2 = reinterpret_cast<const unsigned char *>(buffers[1].data);
                if (ends[0] <= buffers[1].size && (!ends[0] || r2[ends[0] - 1] == '\n') &&
                    lines[1].lines_before(r2, ends[0]) == record * 4) {
                    ends[1] = ends[0]; ++mate_same_offset;
                } else {
                    ends[1] = record ? lines[1].end_of_line(r2, size_t(record * 4)) : 0;
                    ++mate_corrected_offset;
                }
            }
        }
        out[i].first_record = int64_t(previous);
        out[i].records = int64_t(record - previous);
        for (int m = 0; m < mates; ++m) {
            out[i].begin[m] = int64_t(offsets[m]); out[i].end[m] = int64_t(ends[m]);
        }
        // Inspect only partition endpoints here. Full validation and exact PE
        // name comparison happen once, while each CPE copies its FASTQ fields.
        if (record > previous) {
            swbwa_fastq_record_view view[2];
            for (int m = 0; m < mates; ++m) {
                size_t position = offsets[m];
                if (swbwa_fastq_record_next(buffers[m].data, ends[m], &position, &view[m]) != 1)
                    throw std::runtime_error("invalid FASTQ at CPE partition boundary");
            }
            if (mates == 2 && (view[0].name_size != view[1].name_size ||
                std::memcmp(view[0].name, view[1].name, view[0].name_size)))
                throw std::runtime_error("paired FASTQ name token mismatch at partition boundary");
        }
        previous = record; offsets[0] = ends[0]; offsets[1] = ends[1];
    }
}

int swbwa_mpi_gzip_reader::next_static(swbwa_mpi_gzip_batch &out)
{
    if (decoded_done) return 0;
    while (!k) advance(0); // Sampling is Stage 1 work, never open work.
    const uint64_t base_id = multiply(wave, uint64_t(ranks));
    const uint64_t wave_first = multiply(base_id, k);
    const uint64_t wave_end = multiply(add(base_id, uint64_t(ranks)), k);
    for (;;) {
        Stream &a = streams[0], &b = streams[1];
        if (mates == 2 && ((a.eof && b.records > a.records) || (b.eof && a.records > b.records) ||
            (a.eof && b.eof && a.records != b.records)))
            throw std::runtime_error("paired FASTQ record count mismatch");
        const bool eof = a.eof && (mates == 1 || b.eof);
        const uint64_t available = mates == 1 ? a.records : std::min(a.records, b.records);
        if (eof && wave_first >= available) {
            for (int m = 0; m < mates; ++m)
                if (!streams[m].chunks.empty() || streams[m].carry_size)
                    throw std::runtime_error("undelivered MPI FASTQ records at EOF");
            const double wait_begin = MPI_Wtime();
            check(MPI_Barrier(comm));
            eof_wait_s += MPI_Wtime() - wait_begin;
            decoded_done = true; return 0;
        }
        if (available >= wave_end || eof) {
            const double begin = MPI_Wtime();
            const uint64_t id = add(base_id, uint64_t(rank));
            const uint64_t nominal_first = multiply(id, k);
            const uint64_t first = std::min(nominal_first, available);
            const uint64_t count = nominal_first < available ? std::min(k, available - nominal_first) : 0;
            if (id > INT64_MAX || first > uint64_t(INT64_MAX / mates) ||
                count > uint64_t(INT64_MAX / mates) - first)
                throw std::overflow_error("MPI FASTQ n_processed range exceeded");
            fq::OwnedChunk buffers[2];
            SwbwaLineIndex lines[2];
            for (int m = 0; m < mates; ++m) buffers[m] = assemble(streams[m], id, first, count, lines[m]);
            std::unique_ptr<swbwa_fastq_partition_t, decltype(&std::free)> parts(
                count ? static_cast<swbwa_fastq_partition_t *>(
                    std::calloc(workers, sizeof(swbwa_fastq_partition_t))) : NULL, &std::free);
            if (count && !parts) throw std::bad_alloc();
            if (count) partitions(buffers, lines, count, parts.get());
            // Finish boundary planning on every rank before READY. Full record
            // validation is performed by the subsequent CPE formatting pass.
            const double wait_begin = MPI_Wtime();
            if (!dynamic) check(MPI_Barrier(comm));
            ready_wait_s += MPI_Wtime() - wait_begin;
            out.partitions = parts.release(); out.partition_count = count ? workers : 0;
            out.chunk_id = int64_t(id); out.first_record = int64_t(first); out.records = int64_t(count);
            for (int m = 0; m < mates; ++m) {
                out.size[m] = int64_t(buffers[m].size); out.data[m] = buffers[m].release();
            }
            ++wave;
            emit_s += MPI_Wtime() - begin;
            return 1;
        }
        int mate = 0;
        if (mates == 2 && (a.eof || (!b.eof && b.records < a.records))) mate = 1;
        advance(mate);
    }
}

void swbwa_mpi_gzip_reader::release_ready()
{
    if (data_window != MPI_WIN_NULL) {
        check(MPI_Win_unlock_all(data_window));
        check(MPI_Win_free(&data_window));
    }
    ready_payload = fq::OwnedChunk{};
    descriptors.clear(); ready_count = 0;
}

bool swbwa_mpi_gzip_reader::refill()
{
    double begin = MPI_Wtime();
    {
        /* A standard-mode decoder needs the processor thread. Do not enqueue
         * collective codec work ahead of older chunks whose ordered writer
         * can hold another rank's processor at the full output queue. */
        thread_call(pthread_mutex_lock(&pending_mutex));
        while (!pending.empty()) thread_call(pthread_cond_wait(&completed, &pending_mutex));
        thread_call(pthread_mutex_unlock(&pending_mutex));
    }
    drain_wait_s += MPI_Wtime() - begin;
    begin = MPI_Wtime();
    /* Every requester has flushed its Gets before reaching this barrier.
     * Delivered buffers are independent of the immutable window storage. */
    check(MPI_Barrier(comm));
    refill_wait_s += MPI_Wtime() - begin;
    release_ready();
    if (decoded_done) return false;
    const double prepare_begin = MPI_Wtime();
    ready_begin = multiply(wave, uint64_t(ranks));
    uint64_t bytes = 0;
    std::vector<std::unique_ptr<ReadyChunk> > ready;
    std::vector<U64> local_values(size_t(ready_limit) * 9, 0);
    for (unsigned slot = 0; slot < ready_limit; ++slot) {
        std::unique_ptr<ReadyChunk> chunk(new ReadyChunk);
        if (!next_static(chunk->batch)) break;
        const swbwa_mpi_gzip_batch &b = chunk->batch;
        U64 *values = local_values.data() + slot * 9;
        values[0] = U64(b.chunk_id); values[1] = U64(b.first_record);
        values[2] = U64(b.records); values[3] = U64(b.size[0]);
        values[4] = U64(b.size[1]); values[5] = b.partition_count;
        if (b.records) {
            const size_t sizes[3] = {size_t(b.size[0]), size_t(b.size[1]),
                                    b.partition_count * sizeof(*b.partitions)};
            for (unsigned i = 0; i < 3; ++i) {
                values[6 + i] = bytes;
                bytes = add(bytes, sizes[i]);
            }
        }
        ready.push_back(std::move(chunk));
    }
    ready_peak_bytes = std::max(ready_peak_bytes, bytes);
    const size_t slots = size_t(ready_limit) * size_t(ranks);
    std::vector<U64> values(slots * 9, 0);
    check(MPI_Allgather(local_values.data(), int(local_values.size()), MPI_UNSIGNED_LONG_LONG,
        values.data(), int(local_values.size()), MPI_UNSIGNED_LONG_LONG, comm));
    descriptors.resize(slots);
    bool exhausted = false;
    for (size_t index = 0; index < slots; ++index) {
        const size_t source = (index % size_t(ranks)) * ready_limit + index / size_t(ranks);
        ReadyDescriptor &d = descriptors[index];
        std::copy_n(values.data() + source * 9, 9, d.values.begin());
        if (!d.values[2]) { exhausted = true; continue; }
        if (exhausted || d.values[0] != add(ready_begin, index) ||
            d.values[1] != multiply(d.values[0], k) || d.values[5] != workers)
            throw std::runtime_error("noncontiguous MPI gzip ready window");
        ++ready_count;
    }
    if (ready_count) {
        /* SWCH's dynamic-window Get fails with SRC_ADDR_ERR even in an isolated
         * probe. Use one ordinary window with relative byte displacements. */
        if (bytes > SIZE_MAX || bytes > uint64_t(std::numeric_limits<MPI_Aint>::max()))
            throw std::overflow_error("MPI gzip ready payload exceeds address range");
        const double pack_begin = MPI_Wtime();
        ready_payload.data = static_cast<char *>(std::malloc(bytes ? size_t(bytes) : 1));
        if (!ready_payload.data) throw std::bad_alloc();
        ready_payload.size = size_t(bytes);
        for (size_t slot = 0; slot < ready.size(); ++slot) {
            const auto &b = ready[slot]->batch;
            const size_t sizes[3] = {size_t(b.size[0]), size_t(b.size[1]),
                                    b.partition_count * sizeof(*b.partitions)};
            for (unsigned i = 0; i < 3; ++i) if (sizes[i])
                copy_bytes(ready_payload.data + local_values[slot * 9 + 6 + i],
                           static_cast<const unsigned char *>(ready[slot]->buffer(i)), sizes[i]);
            ready[slot].reset();
        }
        pack_s += MPI_Wtime() - pack_begin;
        check(MPI_Win_create(ready_payload.data, MPI_Aint(bytes), 1, MPI_INFO_NULL, comm, &data_window));
        check(MPI_Win_lock_all(0, data_window));
        check(MPI_Win_sync(data_window));
    }
    /* Use RMA for reset too: never race local stores against public counters. */
    if (!rank) {
        const U64 zero = 0;
        check(MPI_Accumulate(&zero, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, 1,
                            MPI_UNSIGNED_LONG_LONG, MPI_REPLACE, ticket_window));
        check(MPI_Win_flush(0, ticket_window));
    }
    check(MPI_Barrier(comm));
    ++ready_windows;
    refill_s += MPI_Wtime() - prepare_begin;
    if (!ready_count && !decoded_done)
        throw std::logic_error("empty ready window before validated gzip EOF");
    return ready_count != 0;
}

int swbwa_mpi_gzip_reader::next_dynamic(swbwa_mpi_gzip_batch &out)
{
    for (;;) {
        if (!ready_count && !refill()) return 0;
        const double begin = MPI_Wtime();
        const U64 one = 1;
        U64 index = 0;
        check(MPI_Fetch_and_op(&one, &index, MPI_UNSIGNED_LONG_LONG, 0, 0,
                              MPI_SUM, ticket_window));
        check(MPI_Win_flush(0, ticket_window));
        claim_s += MPI_Wtime() - begin;
        if (index >= ready_count) {
            if (!refill()) return 0;
            continue;
        }
        const int owner = int(index % uint64_t(ranks));
        const ReadyDescriptor &d = descriptors[size_t(index)];
        const double fetch_begin = MPI_Wtime();
        {
            std::unique_ptr<ReadyChunk> chunk(new ReadyChunk);
            swbwa_mpi_gzip_batch &b = chunk->batch;
            b.chunk_id = int64_t(d.values[0]); b.first_record = int64_t(d.values[1]);
            b.records = int64_t(d.values[2]); b.partition_count = unsigned(d.values[5]);
            for (int m = 0; m < mates; ++m) {
                b.size[m] = int64_t(d.values[3 + m]);
                if (b.size[m] <= 0 || b.size[m] > SWBWA_CPE_FORMAT_BUFFER_BYTES)
                    throw std::runtime_error("invalid MPI gzip ready buffer size");
                b.data[m] = static_cast<char *>(std::malloc(size_t(b.size[m]) + 1));
                if (!b.data[m]) throw std::bad_alloc();
                b.data[m][b.size[m]] = 0;
            }
            b.partitions = static_cast<swbwa_fastq_partition_t *>(
                std::malloc(b.partition_count * sizeof(*b.partitions)));
            if (!b.partitions) throw std::bad_alloc();
            const size_t sizes[3] = {size_t(b.size[0]), size_t(b.size[1]),
                                    b.partition_count * sizeof(*b.partitions)};
            for (unsigned i = 0; i < 3; ++i) {
                if (owner == rank) {
                    if (sizes[i]) copy_bytes(static_cast<char *>(chunk->buffer(i)),
                        reinterpret_cast<const unsigned char *>(ready_payload.data) + d.values[6 + i], sizes[i]);
                    continue;
                }
                size_t offset = 0;
                while (offset < sizes[i]) {
                    const int n = int(std::min<size_t>(INT_MAX, sizes[i] - offset));
                    check(MPI_Get(static_cast<char *>(chunk->buffer(i)) + offset, n, MPI_BYTE,
                        owner, MPI_Aint(d.values[6 + i] + offset), n, MPI_BYTE, data_window));
                    offset += size_t(n);
                }
                fetched_bytes += sizes[i];
            }
            if (owner != rank) check(MPI_Win_flush(owner, data_window));
            out = b; b = swbwa_mpi_gzip_batch{};
        }
        if (owner == rank) {
            ++local_claims;
            local_copy_s += MPI_Wtime() - fetch_begin;
        } else {
            ++remote_claims;
            fetch_s += MPI_Wtime() - fetch_begin;
        }
        return 1;
    }
}

int swbwa_mpi_gzip_reader::next(swbwa_mpi_gzip_batch &out)
{
    thread_check();
    if (done) return 0;
    const int result = dynamic ? next_dynamic(out) : next_static(out);
    if (!result) { done = true; return 0; }
    if (dynamic) {
        thread_call(pthread_mutex_lock(&pending_mutex));
        pending.push_back(out.chunk_id);
        thread_call(pthread_mutex_unlock(&pending_mutex));
    }
    emitted += uint64_t(out.records); emitted_chunks += out.records != 0; empty += out.records == 0;
    for (int m = 0; m < mates; ++m) emitted_bytes[m] += uint64_t(out.size[m]);
    if (debug) std::fprintf(stderr,
        "[MPI gzip rank=%d] chunk=%lld first_record=%lld end_record_exclusive=%lld records=%lld R1_bytes=%lld R2_bytes=%lld\n",
        rank, (long long)out.chunk_id, (long long)out.first_record,
        (long long)(out.first_record + out.records), (long long)out.records,
        (long long)out.size[0], (long long)out.size[1]);
    return 1;
}

extern "C" swbwa_mpi_gzip_reader *swbwa_mpi_gzip_open(const char *path1, const char *path2,
    size_t batch_bytes, unsigned workers, int debug)
{
    swbwa_mpi_gzip_reader *r = NULL;
    last_error[0] = 0;
    try {
        int initialized = 0, finalized = 0;
        check(MPI_Initialized(&initialized));
        if (initialized) check(MPI_Finalized(&finalized));
        if (!initialized || finalized) throw std::runtime_error("MPI gzip requires active MPI");
        thread_check();
        r = new swbwa_mpi_gzip_reader;
        check(MPI_Comm_dup(MPI_COMM_WORLD, &r->comm));
        check(MPI_Comm_set_errhandler(r->comm, MPI_ERRORS_RETURN));
        check(MPI_Comm_rank(r->comm, &r->rank)); check(MPI_Comm_size(r->comm, &r->ranks));
        if (!path1 || !*path1 || (path2 && !*path2) || !batch_bytes ||
            batch_bytes > size_t(SWBWA_CPE_FORMAT_BUFFER_BYTES) || !workers || workers > 384)
            throw std::invalid_argument("invalid MPI gzip reader configuration");
#if !defined(SWBWA_CODEC_RAPID) || !SWBWA_CODEC_RAPID
        throw std::runtime_error("MPI gzip requires rapid backend");
#endif
        r->mates = path2 ? 2 : 1; r->workers = workers; r->debug = debug; r->target = batch_bytes;
        r->ready_limit = r->dynamic ? ready_chunk_limit() : 0;
        r->queue_limit = std::max<uint64_t>(512U << 20, multiply(batch_bytes, 3));
        std::array<U64, 6> config = {{U64(r->mates), batch_bytes, workers, U64(debug), U64(r->dynamic), r->ready_limit}};
        std::vector<std::array<U64, 6> > all(size_t(r->ranks), std::array<U64, 6>{{}});
        check(MPI_Allgather(config.data(), 6, MPI_UNSIGNED_LONG_LONG, all.data(), 6, MPI_UNSIGNED_LONG_LONG, r->comm));
        for (const auto &c : all) if (c != config) throw std::runtime_error("MPI gzip reader configuration mismatch");
        swbwa_mpi_decoder_config decoder_config;
        swbwa_mpi_decoder_config_init(&decoder_config); decoder_config.workers = workers;
        for (int m = 0; m < r->mates; ++m) {
            char error[512] = "";
            if (swbwa_mpi_decoder_open(r->comm, m ? path2 : path1, &decoder_config,
                                       &r->streams[m].decoder, error, sizeof(error)))
                throw std::runtime_error(error);
        }
        if (r->dynamic) {
            check(MPI_Win_create(&r->ticket, r->rank == 0 ? sizeof(r->ticket) : 0,
                                 sizeof(r->ticket), MPI_INFO_NULL, r->comm, &r->ticket_window));
            check(MPI_Win_lock_all(0, r->ticket_window));
        }
        return r;
    } catch (const std::exception &e) { abort_reader(r, e.what()); }
    catch (...) { abort_reader(r, "unknown MPI gzip open error"); }
    delete r; return NULL;
}

extern "C" int swbwa_mpi_gzip_next(swbwa_mpi_gzip_reader *r, swbwa_mpi_gzip_batch *out)
{
    if (out) *out = swbwa_mpi_gzip_batch{};
    try {
        if (!r || !out) throw std::invalid_argument("null MPI gzip reader/batch");
        const double begin = MPI_Wtime();
        const int status = r->next(*out);
        r->next_s += MPI_Wtime() - begin;
        return status;
    } catch (const std::exception &e) { return abort_reader(r, e.what()); }
    catch (...) { return abort_reader(r, "unknown MPI gzip next error"); }
}

extern "C" const char *swbwa_mpi_gzip_error(const swbwa_mpi_gzip_reader *r)
{
    return r ? r->error : last_error;
}

extern "C" int swbwa_mpi_gzip_complete(swbwa_mpi_gzip_reader *r, int64_t id)
{
    try {
        if (!r) throw std::invalid_argument("null MPI gzip completion reader");
        if (!r->dynamic) return 0;
        thread_call(pthread_mutex_lock(&r->pending_mutex));
        if (r->pending.empty() || r->pending.front() != id) {
            thread_call(pthread_mutex_unlock(&r->pending_mutex));
            throw std::logic_error("MPI gzip completion out of order or repeated");
        }
        r->pending.pop_front();
        thread_call(pthread_cond_signal(&r->completed));
        thread_call(pthread_mutex_unlock(&r->pending_mutex));
        return 0;
    } catch (const std::exception &e) { return abort_reader(r, e.what()); }
    catch (...) { return abort_reader(r, "unknown MPI gzip completion error"); }
}

extern "C" int64_t swbwa_mpi_gzip_chunk_count(const swbwa_mpi_gzip_reader *r)
{
    if (!r || !r->done || !r->k) return -1;
    const uint64_t records = r->streams[0].records;
    const uint64_t chunks = records / r->k + (records % r->k != 0);
    return chunks <= uint64_t(INT64_MAX) ? int64_t(chunks) : -1;
}

extern "C" int swbwa_mpi_gzip_report(swbwa_mpi_gzip_reader *r, FILE *file)
{
    if (!r || !file) return abort_reader(r, "null MPI gzip report argument");
    const int result = std::fprintf(file,
        "[MPI gzip rank=%d] records=%llu chunks=%llu empty_READY=%llu records_per_chunk=%llu "
        "R1_decoded_bytes=%llu R2_decoded_bytes=%llu next_s=%.6f decode_s=%.6f frame_carry_s=%.6f "
        "exchange_s=%.6f emit_layout_s=%.6f validated_EOF=%d queue_limit_per_mate_bytes=%llu\n",
        r->rank, U64(r->emitted), U64(r->emitted_chunks), U64(r->empty), U64(r->k),
        U64(r->emitted_bytes[0]), U64(r->emitted_bytes[1]), r->next_s, r->decode_s,
        r->frame_s, r->exchange_s, r->emit_s, int(r->done), U64(r->queue_limit));
    if (result < 0) return abort_reader(r, "MPI gzip report write failed");
    if (std::fprintf(file,
        "[MPI gzip rank=%d scheduler] mode=%s ready_chunks_per_rank=%u windows=%llu "
        "local_claims=%llu remote_claims=%llu fetched_bytes=%llu ready_peak_bytes=%llu "
        "refill_s=%.6f refill_wait_s=%.6f drain_wait_s=%.6f claim_rma_s=%.6f fetch_s=%.6f "
        "pack_s=%.6f local_copy_s=%.6f\n",
        r->rank, r->dynamic ? "dynamic" : "static", r->ready_limit, U64(r->ready_windows),
        U64(r->local_claims), U64(r->remote_claims), U64(r->fetched_bytes), U64(r->ready_peak_bytes),
        r->refill_s, r->refill_wait_s, r->drain_wait_s, r->claim_s, r->fetch_s,
        r->pack_s, r->local_copy_s) < 0)
        return abort_reader(r, "MPI gzip scheduler report write failed");
    if (std::fprintf(file,
        "[MPI gzip rank=%d reader_local] frame_local_s=%.6f exchange_local_s=%.6f "
        "emit_local_s=%.6f adapter_local_s=%.6f frame_mpi_s=%.6f exchange_mpi_s=%.6f "
        "ready_wait_s=%.6f eof_wait_s=%.6f local_payload_bytes=%llu remote_payload_bytes=%llu "
        "assembled_copy_bytes=%llu mate_same_offset=%llu mate_corrected_offset=%llu "
        "record_validation=CPE_format timers=rank_local_wall_seconds\n",
        r->rank, r->frame_s - r->frame_mpi_s, r->exchange_s - r->exchange_mpi_s,
        r->emit_s - r->ready_wait_s,
        r->frame_s + r->exchange_s + r->emit_s - r->frame_mpi_s - r->exchange_mpi_s - r->ready_wait_s,
        r->frame_mpi_s, r->exchange_mpi_s, r->ready_wait_s, r->eof_wait_s,
        U64(r->local_payload_bytes), U64(r->remote_payload_bytes), U64(r->assembled_copy_bytes),
        U64(r->mate_same_offset), U64(r->mate_corrected_offset)) < 0)
        return abort_reader(r, "MPI gzip report write failed");
    for (int m = 0; m < r->mates; ++m) {
        swbwa_mpi_decoder_metrics d = {};
        if (swbwa_mpi_decoder_get_metrics(r->streams[m].decoder, &d))
            return abort_reader(r, "MPI gzip decoder metrics unavailable");
        if (std::fprintf(file,
            "[MPI gzip rank=%d mate=%d] epochs=%llu decoded_bytes=%llu source_read_bytes=%llu "
            "io_s=%.6f discover_run_s=%.6f discover_wait_s=%.6f resolve_run_s=%.6f resolve_wait_s=%.6f "
            "queue_peak_accounted_bytes=%llu boundary_records=%llu\n", r->rank, m + 1, U64(d.epochs), U64(d.decoded_bytes),
            U64(d.source_read_bytes), d.io_s, d.discover_run_s, d.discover_wait_s, d.resolve_run_s,
            d.resolve_wait_s, U64(r->streams[m].peak_queued), U64(r->streams[m].boundary_records)) < 0)
            return abort_reader(r, "MPI gzip report write failed");
        const double next_io = d.io_s - d.open_io_s;
        const double accounted = next_io + d.prepare_s + d.cleanup_s + d.discover_wall_s +
            d.resolve_wall_s + d.metadata_local_s + d.metadata_mpi_s + d.window_local_s +
            d.window_mpi_s + d.commit_local_s + d.commit_mpi_s + d.commit_barrier_s + d.eof_barrier_s;
        if (std::fprintf(file,
            "[MPI gzip rank=%d mate=%d decoder] open_s=%.6f open_io_s=%.6f next_s=%.6f "
            "next_io_s=%.6f prepare_s=%.6f cleanup_s=%.6f other_s=%.6f\n"
            "[MPI gzip rank=%d mate=%d decoder] metadata_local_s=%.6f metadata_mpi_s=%.6f "
            "window_local_s=%.6f window_mpi_s=%.6f commit_local_s=%.6f commit_mpi_s=%.6f "
            "commit_barrier_s=%.6f eof_barrier_s=%.6f\n"
            "[MPI gzip rank=%d mate=%d decoder] discover_wall_s=%.6f resolve_wall_s=%.6f "
            "capacity_retries=%llu retry_s=%.6f retry_jobs=%llu\n"
            "[MPI gzip rank=%d mate=%d decoder] discover_jobs=%llu accepted_jobs=%llu "
            "speculative_accepted=%llu discarded_jobs=%llu discarded_output_bytes=%llu "
            "discover_failed_jobs=%llu discover_output_full_jobs=%llu peak_buffer_bytes=%llu\n"
            "[MPI gzip rank=%d mate=%d decoder] timers=rank_local_wall_seconds "
            "mpi_includes_peer_wait=1 run_wait_included_in_cpe_wall=1 retry_included_in_next=1 "
            "job_accounting=discover:accepted+discarded+retry\n",
            r->rank, m + 1, d.open_s, d.open_io_s, d.next_s, next_io, d.prepare_s, d.cleanup_s,
            std::max(0.0, d.next_s - accounted),
            r->rank, m + 1, d.metadata_local_s, d.metadata_mpi_s, d.window_local_s, d.window_mpi_s,
            d.commit_local_s, d.commit_mpi_s, d.commit_barrier_s, d.eof_barrier_s,
            r->rank, m + 1, d.discover_wall_s, d.resolve_wall_s, U64(d.capacity_retries),
            d.retry_s, U64(d.retry_jobs),
            r->rank, m + 1, U64(d.discover_jobs), U64(d.accepted_jobs), U64(d.speculative_accepted),
            U64(d.discarded_jobs), U64(d.discarded_output_bytes), U64(d.discover_failed_jobs),
            U64(d.discover_output_full_jobs), U64(d.peak_buffer_bytes), r->rank, m + 1) < 0)
            return abort_reader(r, "MPI gzip report write failed");
    }
    return 0;
}

extern "C" void swbwa_mpi_gzip_close(swbwa_mpi_gzip_reader *r)
{
    if (!r) return;
    try {
        thread_check();
        if (!r->done) throw std::runtime_error("MPI gzip close requires globally validated EOF from next");
        if (r->dynamic) {
            check(MPI_Barrier(r->comm));
            r->release_ready();
            check(MPI_Win_unlock_all(r->ticket_window));
            check(MPI_Win_free(&r->ticket_window));
        }
        for (int m = 0; m < r->mates; ++m) {
            if (!r->streams[m].eof || r->streams[m].carry_size || !r->streams[m].chunks.empty())
                throw std::runtime_error("MPI gzip close before mate CRC/framing completion");
            if (swbwa_mpi_decoder_close(r->streams[m].decoder) != 0)
                throw std::runtime_error(swbwa_mpi_decoder_error(r->streams[m].decoder));
        }
        check(MPI_Comm_free(&r->comm)); delete r;
    } catch (const std::exception &e) { abort_reader(r, e.what()); }
    catch (...) { abort_reader(r, "unknown MPI gzip close error"); }
}
