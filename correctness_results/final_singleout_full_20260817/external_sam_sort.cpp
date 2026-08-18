#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Record {
    std::uint64_t id;
    std::uint64_t order;
    std::string line;
};

static std::uint64_t parse_id(const std::string &line)
{
    const std::size_t qname_end = line.find_first_of("\t ");
    const std::size_t end = qname_end == std::string::npos
                          ? line.size() : qname_end;
    const std::size_t dot = end == 0 ? std::string::npos
                                     : line.rfind('.', end - 1);
    std::uint64_t value = 0;

    if (dot == std::string::npos || dot + 1 >= end)
        throw std::runtime_error("invalid SAM QNAME: " + line);

    for (std::size_t i = dot + 1; i < end; ++i) {
        const unsigned digit = static_cast<unsigned>(line[i] - '0');

        if (digit > 9)
            throw std::runtime_error("non-numeric SAM QNAME suffix: " + line);
        if (value > (std::numeric_limits<std::uint64_t>::max() - digit) / 10)
            throw std::runtime_error("SAM QNAME suffix overflow: " + line);
        value = value * 10 + digit;
    }
    return value;
}

static std::string run_path(const std::string &directory, std::size_t index)
{
    std::ostringstream path;

    path << directory << "/swbwa_sort_run_"
         << std::setw(6) << std::setfill('0') << index << ".tmp";
    return path.str();
}

static void write_run(std::vector<Record> *records, const std::string &path)
{
    std::sort(records->begin(), records->end(),
              [](const Record &a, const Record &b) {
                  return a.id < b.id || (a.id == b.id && a.order < b.order);
              });

    std::ofstream output(path.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
        throw std::runtime_error("cannot create temporary run: " + path);

    for (std::size_t i = 0; i < records->size(); ++i)
        output << (*records)[i].line << '\n';
    output.close();
    if (!output)
        throw std::runtime_error("failed to write temporary run: " + path);
}

struct RunReader {
    explicit RunReader(const std::string &path) : input(path.c_str(), std::ios::binary) {}

    std::ifstream input;
    std::string line;
    std::uint64_t id;
};

struct HeapEntry {
    std::uint64_t id;
    std::size_t run;
};

struct HeapGreater {
    bool operator()(const HeapEntry &a, const HeapEntry &b) const
    {
        return a.id > b.id || (a.id == b.id && a.run > b.run);
    }
};

static bool read_next(RunReader *reader)
{
    if (!std::getline(reader->input, reader->line)) {
        if (reader->input.eof())
            return false;
        throw std::runtime_error("failed to read a temporary run");
    }
    if (!reader->line.empty() && reader->line.back() == '\r')
        reader->line.pop_back();
    reader->id = parse_id(reader->line);
    return true;
}

static void merge_runs(const std::vector<std::string> &paths,
                       const std::string &output_path)
{
    std::vector<std::unique_ptr<RunReader> > readers;
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapGreater> heap;

    readers.reserve(paths.size());
    for (std::size_t i = 0; i < paths.size(); ++i) {
        readers.push_back(std::unique_ptr<RunReader>(new RunReader(paths[i])));
        if (!readers.back()->input)
            throw std::runtime_error("cannot open temporary run: " + paths[i]);
        if (read_next(readers.back().get()))
            heap.push(HeapEntry{readers.back()->id, i});
    }

    std::ofstream output(output_path.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
        throw std::runtime_error("cannot create output: " + output_path);

    std::uint64_t written = 0;
    while (!heap.empty()) {
        const HeapEntry entry = heap.top();
        RunReader *reader = readers[entry.run].get();

        heap.pop();
        output << reader->line << '\n';
        if (++written % 1000000 == 0)
            std::cerr << "Merged " << written << " SAM records\n";
        if (read_next(reader))
            heap.push(HeapEntry{reader->id, entry.run});
    }
    output.close();
    if (!output)
        throw std::runtime_error("failed to write output: " + output_path);
}

int main(int argc, char **argv)
{
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << argv[0]
                  << " input.sam output.sam temporary-directory [chunk-MiB]\n";
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const std::string temp_directory = argv[3];
    const std::size_t chunk_mib = argc == 5
                                ? static_cast<std::size_t>(std::strtoull(argv[4], NULL, 10))
                                : 128;
    const std::size_t chunk_bytes = chunk_mib * 1024 * 1024;
    std::vector<Record> records;
    std::vector<std::string> paths;
    std::uint64_t order = 0;
    std::size_t buffered_bytes = 0;

    if (input_path == output_path || chunk_mib == 0 ||
        chunk_mib > std::numeric_limits<std::size_t>::max() / (1024 * 1024)) {
        std::cerr << "Invalid input/output path or chunk size\n";
        return 1;
    }

    try {
        std::ifstream input(input_path.c_str(), std::ios::binary);
        std::string line;

        if (!input)
            throw std::runtime_error("cannot open input: " + input_path);

        while (std::getline(input, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            buffered_bytes += line.size() + sizeof(Record) + 1;
            records.push_back(Record{parse_id(line), order++, std::move(line)});

            if (buffered_bytes >= chunk_bytes) {
                const std::string path = run_path(temp_directory, paths.size());
                write_run(&records, path);
                paths.push_back(path);
                std::cerr << "Wrote run " << paths.size()
                          << ", total records " << order << '\n';
                records.clear();
                buffered_bytes = 0;
            }
        }
        if (!input.eof()) {
            std::ostringstream message;
            message << "failed to read input: state=" << input.rdstate()
                    << " errno=" << errno << " (" << std::strerror(errno) << ')';
            throw std::runtime_error(message.str());
        }
        if (!records.empty()) {
            const std::string path = run_path(temp_directory, paths.size());
            write_run(&records, path);
            paths.push_back(path);
            std::cerr << "Wrote run " << paths.size()
                      << ", total records " << order << '\n';
        }
        if (paths.empty()) {
            std::ofstream empty(output_path.c_str(), std::ios::binary | std::ios::trunc);
            if (!empty)
                throw std::runtime_error("cannot create empty output: " + output_path);
        } else {
            merge_runs(paths, output_path);
        }

        for (std::size_t i = 0; i < paths.size(); ++i)
            std::remove(paths[i].c_str());
        std::cout << "Sorted " << order << " SAM records into "
                  << output_path << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 2;
    }
}
