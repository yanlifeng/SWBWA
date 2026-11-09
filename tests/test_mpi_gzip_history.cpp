#include "mpi_history_tree.h"
#include <cstdio>
#include <random>

using namespace swbwa_mpi_gzip;

static void require(bool value, const char *text) {
    if (!value) throw std::runtime_error(text);
}

static void run(MPI_Comm comm) {
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    for (unsigned fixture = 0; fixture < 4; ++fixture) {
        std::mt19937 rng(749 + fixture);
        Window initial(fixture == 0 ? 0 : fixture == 1 ? 17 : HISTORY);
        for (auto &v : initial) v = rng() & 255;
        Window current = initial;
        std::vector<HistorySummary> summaries(ranks);
        std::vector<Window> incoming(ranks);
        for (int q = 0; q < ranks; ++q) {
            incoming[q] = current;
            std::vector<std::vector<uint16_t>> jobs;
            for (unsigned j = 0; j < (fixture == 3 ? 0 : 3); ++j) {
                const size_t n = q % 4 == 2 ? 0 :
                    fixture == 1 ? 9 + rng() % 43 : rng() % (HISTORY * 2);
                jobs.emplace_back(n);
                for (auto &symbol : jobs.back()) {
                    symbol = current.empty() || rng() % 3 == 0 ? rng() % 256 :
                        65536U - current.size() + rng() % current.size();
                }
                current = advance(current, jobs.back().data(), jobs.back().size());
            }
            std::vector<SymbolSpan> spans;
            for (const auto &job : jobs) spans.push_back({job.data(), job.size()});
            summarize_history(spans, summaries[q]);
            require(resolve_history(incoming[q], summaries[q]) == current, "summary oracle");
        }
        HistorySummary identity, total;
        summarize_history({}, identity);
        total = identity;
        for (int q = 0; q < ranks; ++q) {
            HistorySummary out, left_identity, right_identity;
            compose_history(identity, summaries[q], left_identity);
            compose_history(summaries[q], identity, right_identity);
            require(left_identity.symbols == summaries[q].symbols &&
                    right_identity.symbols == summaries[q].symbols &&
                    left_identity.bytes == summaries[q].bytes &&
                    right_identity.bytes == summaries[q].bytes, "identity");
            compose_history(total, summaries[q], out);
            total = out;
        }
        require(resolve_history(initial, total) == current, "composed oracle");
        if (ranks >= 3) {
            HistorySummary ab, bc, abc, abc2;
            compose_history(summaries[0], summaries[1], ab);
            compose_history(summaries[1], summaries[2], bc);
            compose_history(ab, summaries[2], abc);
            compose_history(summaries[0], bc, abc2);
            require(abc.symbols == abc2.symbols && abc.bytes == abc2.bytes, "associativity");
        }
        for (unsigned anchor = 0; anchor < unsigned(ranks); ++anchor) {
            swbwa_mpi_decoder_metrics metrics{};
            HistoryTree tree(comm, rank, ranks, anchor, metrics);
            const unsigned logical = (rank + ranks - anchor) % ranks;
            tree.local_summary() = summaries[logical];
            Window actual, final;
            tree.propagate(initial, actual, final);
            require(actual == incoming[logical], "tree incoming differs");
            require(final == current, "tree final differs");
            uint64_t counts[3] = {metrics.history_tree_compose_calls,
                metrics.history_tree_apply_calls, metrics.history_tree_messages};
            uint64_t sums[3] = {};
            MPI_Allreduce(counts, sums, 3, MPI_UINT64_T, MPI_SUM, comm);
            require(sums[0] == unsigned(ranks - 1) && sums[1] == unsigned(ranks) &&
                    sums[2] == unsigned(2 * (ranks - 1)), "tree operation count");
        }
    }
    HistorySummary identity, invalid, out;
    summarize_history({}, identity);
    invalid = identity;
    invalid.symbols[0] = 256;
    bool rejected = false;
    try { compose_history(identity, invalid, out); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected, "invalid marker accepted");
    rejected = false;
    try { compose_history(identity, identity, identity); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected, "alias accepted");
    invalid = identity; invalid.bytes = UINT64_MAX; identity.bytes = 1;
    rejected = false;
    try { compose_history(identity, invalid, out); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected, "overflow accepted");
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try {
        MPI_Comm comm;
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
        run(comm);
        MPI_Comm_free(&comm);
        MPI_Comm_split(MPI_COMM_WORLD, 0, ranks - rank, &comm);
        run(comm);
        MPI_Comm_free(&comm);
        if (!rank) std::printf("HISTORY_TREE_PASS ranks=%d\n", ranks);
    } catch (const std::exception &e) {
        std::fprintf(stderr, "rank %d: %s\n", rank, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
    return 0;
}
