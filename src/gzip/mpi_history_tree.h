#ifndef SWBWA_MPI_HISTORY_TREE_H
#define SWBWA_MPI_HISTORY_TREE_H

#include "mpi_decoder.h"
#include "mpi_decoder_protocol.h"
#include <string>

namespace swbwa_mpi_gzip {

// Each interval is owned by its logical leftmost rank. Only child owners
// exchange messages; recursion visits this rank's path, not the whole tree.
// The caller must abort its parent communicator on exceptions, as for decoder
// failures. Epochs must not overlap on this private communicator.
class HistoryTree {
public:
    HistoryTree(MPI_Comm comm, int rank, int ranks, unsigned anchor,
                swbwa_mpi_decoder_metrics &metrics)
        : comm_(comm), ranks_(unsigned(ranks)), anchor_(anchor),
          logical_(logical_rank(rank, ranks, anchor)), metrics_(metrics),
          buffers_(history_tree_depth(ranks_) + 2) {}

    HistorySummary &local_summary() { return buffers_[0]; }

    void propagate(const Window &initial, Window &incoming, Window &final) {
        up(0, ranks_, 0);
        if (!logical_) {
            incoming = initial;
            final = apply(initial, buffers_[0]);
        }
        down(0, ranks_, 0, incoming);
        if (ranks_ == 1) return;
        const double begin = MPI_Wtime();
        unsigned size = unsigned(final.size());
        check(MPI_Bcast(&size, 1, MPI_UNSIGNED, physical(0), comm_));
        if (size > HISTORY) throw std::runtime_error("invalid tree final dictionary size");
        final.resize(size);
        check(MPI_Bcast(final.data(), int(size), MPI_BYTE, physical(0), comm_));
        const double seconds = MPI_Wtime() - begin;
        metrics_.history_final_sync_s += seconds;
        metrics_.window_mpi_s += seconds;
    }

private:
    enum { SUMMARY_TAG = 27101, WINDOW_TAG = 27102 };
    MPI_Comm comm_;
    unsigned ranks_, anchor_, logical_;
    swbwa_mpi_decoder_metrics &metrics_;
    // Aggregate, receive scratch, and one saved left summary per tree level.
    std::vector<HistorySummary> buffers_;

    static unsigned logical_rank(int rank, int ranks, unsigned anchor) {
        if (ranks <= 0 || rank < 0 || rank >= ranks || anchor >= unsigned(ranks))
            throw std::runtime_error("invalid history tree rank");
        return (unsigned(rank) + unsigned(ranks) - anchor) % unsigned(ranks);
    }

    int physical(unsigned logical) const { return int((anchor_ + logical) % ranks_); }

    static void check(int status) {
        if (status == MPI_SUCCESS) return;
        char text[MPI_MAX_ERROR_STRING];
        int size = 0;
        MPI_Error_string(status, text, &size);
        throw std::runtime_error(std::string(text, size));
    }

    void send(const void *data, int size, unsigned destination, int tag) {
        const double begin = MPI_Wtime();
        check(MPI_Send(data, size, MPI_BYTE, physical(destination), tag, comm_));
        const double seconds = MPI_Wtime() - begin;
        metrics_.history_tree_mpi_s += seconds;
        metrics_.window_mpi_s += seconds;
        metrics_.history_tree_sent_bytes += unsigned(size);
        ++metrics_.history_tree_messages;
    }

    int receive(void *data, int capacity, unsigned source, int tag) {
        const double begin = MPI_Wtime();
        MPI_Status status;
        check(MPI_Recv(data, capacity, MPI_BYTE, physical(source), tag, comm_, &status));
        int size = 0;
        check(MPI_Get_count(&status, MPI_BYTE, &size));
        const double seconds = MPI_Wtime() - begin;
        metrics_.history_tree_mpi_s += seconds;
        metrics_.window_mpi_s += seconds;
        if (size < 0 || size > capacity) throw std::runtime_error("invalid history message size");
        return size;
    }

    Window apply(const Window &input, const HistorySummary &summary) {
        const double begin = MPI_Wtime();
        Window output = resolve_history(input, summary);
        metrics_.history_tree_apply_s += MPI_Wtime() - begin;
        ++metrics_.history_tree_apply_calls;
        return output;
    }

    void up(unsigned first, unsigned last, unsigned depth) {
        if (last - first == 1) return;
        const unsigned middle = first + (last - first) / 2;
        if (logical_ < middle) {
            up(first, middle, depth + 1);
            if (logical_ != first) return;
            buffers_[depth + 2] = buffers_[0];
            const int size = receive(&buffers_[1], sizeof(HistorySummary), middle, SUMMARY_TAG);
            if (size != int(sizeof(HistorySummary))) throw std::runtime_error("short history summary");
            const double begin = MPI_Wtime();
            compose_history(buffers_[depth + 2], buffers_[1], buffers_[0]);
            metrics_.history_tree_compose_s += MPI_Wtime() - begin;
            ++metrics_.history_tree_compose_calls;
        } else {
            up(middle, last, depth + 1);
            if (logical_ == middle) send(&buffers_[0], sizeof(HistorySummary), first, SUMMARY_TAG);
        }
    }

    void down(unsigned first, unsigned last, unsigned depth, Window &incoming) {
        if (last - first == 1) return;
        const unsigned middle = first + (last - first) / 2;
        if (logical_ == first) {
            Window right = apply(incoming, buffers_[depth + 2]);
            send(right.data(), int(right.size()), middle, WINDOW_TAG);
        } else if (logical_ == middle) {
            incoming.resize(HISTORY);
            incoming.resize(receive(incoming.data(), int(HISTORY), first, WINDOW_TAG));
        }
        if (logical_ < middle) down(first, middle, depth + 1, incoming);
        else down(middle, last, depth + 1, incoming);
    }
};
} // namespace swbwa_mpi_gzip
#endif
