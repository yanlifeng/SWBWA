#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct SamRecord {
    std::uint64_t id;
    std::size_t original_order;
    std::string line;
};

static std::uint64_t parse_name_id(const std::string &line,
                                   std::size_t line_no)
{
    // QNAME 是 SAM 的第一列，以 tab 或空格结束。
    const std::size_t qname_end = line.find_first_of("\t ");
    const std::size_t end =
        (qname_end == std::string::npos) ? line.size() : qname_end;

    if (end == 0) {
        throw std::runtime_error(
            "empty QNAME at line " + std::to_string(line_no));
    }

    // 在 QNAME 范围内寻找最后一个 '.'。
    const std::size_t dot = line.rfind('.', end - 1);

    if (dot == std::string::npos || dot + 1 >= end) {
        throw std::runtime_error(
            "cannot find numeric suffix after '.' at line " +
            std::to_string(line_no) +
            ": " + line);
    }

    std::uint64_t value = 0;
    std::size_t p = dot + 1;
    bool has_digit = false;

    while (p < end && line[p] >= '0' && line[p] <= '9') {
        has_digit = true;

        const unsigned digit =
            static_cast<unsigned>(line[p] - '0');

        if (value >
            (std::numeric_limits<std::uint64_t>::max() - digit) / 10) {
            throw std::runtime_error(
                "numeric ID overflows uint64 at line " +
                std::to_string(line_no));
        }

        value = value * 10 + digit;
        ++p;
    }

    if (!has_digit) {
        throw std::runtime_error(
            "numeric suffix is empty at line " +
            std::to_string(line_no) +
            ": " + line);
    }

    return value;
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr
            << "Usage: " << argv[0]
            << " input.sam output.sam\n";
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];

    if (input_path == output_path) {
        std::cerr
            << "Error: input and output files must be different.\n";
        return 1;
    }

    try {
        std::ifstream input(input_path, std::ios::binary);

        if (!input) {
            throw std::runtime_error(
                "cannot open input file: " + input_path);
        }

        std::vector<SamRecord> records;
        std::string line;
        std::size_t line_no = 0;

        while (std::getline(input, line)) {
            ++line_no;

            // 兼容 Windows CRLF 文件。
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }

            const std::uint64_t id =
                parse_name_id(line, line_no);

            records.push_back(SamRecord{
                id,
                records.size(),
                std::move(line)
            });

            if (line_no % 1000000 == 0) {
                std::cerr
                    << "Loaded " << line_no
                    << " SAM records\n";
            }
        }

        if (!input.eof()) {
            throw std::runtime_error(
                "error while reading input file");
        }

        std::cerr
            << "Finished loading " << records.size()
            << " SAM records\n";

        std::cerr << "Sorting records by ID...\n";

        std::stable_sort(
            records.begin(),
            records.end(),
            [](const SamRecord &a, const SamRecord &b) {
                return a.id < b.id;
            });

        std::cerr << "Sorting finished\n";

        std::ofstream output(
            output_path,
            std::ios::binary | std::ios::trunc);

        if (!output) {
            throw std::runtime_error(
                "cannot open output file: " + output_path);
        }

        for (std::size_t i = 0; i < records.size(); ++i) {
            output << records[i].line << '\n';

            if (!output) {
                throw std::runtime_error(
                    "error while writing output file");
            }

            if ((i + 1) % 1000000 == 0) {
                std::cerr
                    << "Written " << i + 1
                    << " SAM records\n";
            }
        }

        output.close();

        if (!output) {
            throw std::runtime_error(
                "failed to close output file correctly");
        }

        std::cout
            << "Sorted " << records.size()
            << " SAM records\n"
            << "Output: " << output_path << "\n";

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}
