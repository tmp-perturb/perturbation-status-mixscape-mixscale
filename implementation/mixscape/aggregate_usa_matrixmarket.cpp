// Stream an adapter MatrixMarket USA file into gene-level S+A MatrixMarket.

#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include <zlib.h>

namespace {

struct Header {
    long long features{};
    long long cells{};
    long long entries{};
};

struct Parser {
    const char* cursor;
    const char* end;

    bool integer(long long& value) {
        while (cursor < end && (*cursor == ' ' || *cursor == '\t')) {
            ++cursor;
        }
        if (cursor == end) {
            return false;
        }
        auto result = std::from_chars(cursor, end, value);
        if (result.ec != std::errc()) {
            return false;
        }
        cursor = result.ptr;
        return true;
    }
};

Header read_header(gzFile input) {
    char line[512];
    if (gzgets(input, line, sizeof(line)) == nullptr ||
        std::strncmp(line, "%%MatrixMarket matrix coordinate integer general", 48) != 0) {
        throw std::runtime_error("Unexpected MatrixMarket header");
    }
    if (gzgets(input, line, sizeof(line)) == nullptr || line[0] != '%') {
        throw std::runtime_error("Missing MatrixMarket comment line");
    }
    if (gzgets(input, line, sizeof(line)) == nullptr) {
        throw std::runtime_error("Missing MatrixMarket dimensions");
    }
    Parser parser{line, line + std::strlen(line)};
    Header header;
    if (!parser.integer(header.features) || !parser.integer(header.cells) ||
        !parser.integer(header.entries)) {
        throw std::runtime_error("Invalid MatrixMarket dimensions");
    }
    return header;
}

template <typename Callback>
void scan_entries(const std::string& path, long long n_genes, Callback callback) {
    gzFile input = gzopen(path.c_str(), "rb");
    if (input == nullptr) {
        throw std::runtime_error("Cannot open input MatrixMarket: " + path);
    }
    Header header = read_header(input);
    if (header.features != 3 * n_genes) {
        gzclose(input);
        throw std::runtime_error("Input feature count is not 3 * n_genes");
    }
    char line[256];
    long long row, column, value;
    while (gzgets(input, line, sizeof(line)) != nullptr) {
        Parser parser{line, line + std::strlen(line)};
        if (!parser.integer(row) || !parser.integer(column) || !parser.integer(value)) {
            gzclose(input);
            throw std::runtime_error("Invalid MatrixMarket entry");
        }
        callback(row, column, value);
    }
    int error_number = 0;
    const char* error_message = gzerror(input, &error_number);
    if (error_number != Z_OK && error_number != Z_STREAM_END) {
        std::string message = error_message == nullptr ? "gzip read error" : error_message;
        gzclose(input);
        throw std::runtime_error(message);
    }
    gzclose(input);
}

template <typename Callback>
long long aggregate_cells(const std::string& path, long long n_genes, Callback callback) {
    std::vector<long long> values(static_cast<size_t>(n_genes), 0);
    std::vector<long long> touched;
    touched.reserve(12000);
    long long current_column = 0;
    long long aggregated_nnz = 0;
    long long input_nnz = 0;

    auto finish_cell = [&]() {
        for (long long gene : touched) {
            if (values[static_cast<size_t>(gene)] != 0) {
                callback(gene + 1, current_column, values[static_cast<size_t>(gene)]);
                ++aggregated_nnz;
            }
            values[static_cast<size_t>(gene)] = 0;
        }
        touched.clear();
    };

    scan_entries(path, n_genes, [&](long long row, long long column, long long value) {
        ++input_nnz;
        if (input_nnz % 50000000 == 0) {
            std::cerr << "Processed input entries: " << input_nnz << "\n";
        }
        if (row < 1 || row > 3 * n_genes || column < 1 || value < 0) {
            throw std::runtime_error("MatrixMarket entry is outside the expected count range");
        }
        if (current_column == 0) {
            current_column = column;
        }
        if (column < current_column) {
            throw std::runtime_error("MatrixMarket columns are not in nondecreasing order");
        }
        if (column != current_column) {
            finish_cell();
            current_column = column;
        }
        if (row <= n_genes || row > 2 * n_genes) {
            long long gene = (row - 1) % n_genes;
            if (values[static_cast<size_t>(gene)] == 0) {
                touched.push_back(gene);
            }
            values[static_cast<size_t>(gene)] += value;
        }
    });
    if (current_column != 0) {
        finish_cell();
    }
    return aggregated_nnz;
}

class BufferedGzipWriter {
public:
    explicit BufferedGzipWriter(const std::string& path) : output_(gzopen(path.c_str(), "wb1")) {
        if (output_ == nullptr) {
            throw std::runtime_error("Cannot open output MatrixMarket: " + path);
        }
        buffer_.reserve(1 << 20);
    }

    ~BufferedGzipWriter() {
        if (output_ != nullptr) {
            gzclose(output_);
        }
    }

    void text(const std::string& value) {
        buffer_ += value;
        if (buffer_.size() >= (1 << 20)) {
            flush();
        }
    }

    void entry(long long row, long long column, long long value) {
        char line[128];
        int length = std::snprintf(line, sizeof(line), "%lld %lld %lld\n", row, column, value);
        text(std::string(line, static_cast<size_t>(length)));
    }

    void flush() {
        if (!buffer_.empty()) {
            int written = gzwrite(output_, buffer_.data(), static_cast<unsigned>(buffer_.size()));
            if (written != static_cast<int>(buffer_.size())) {
                throw std::runtime_error("gzip write error");
            }
            buffer_.clear();
        }
    }

private:
    gzFile output_;
    std::string buffer_;
};

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: aggregate_usa_matrixmarket <input.mtx.gz> <output.mtx.gz> <n_genes> <n_cells>\n";
        return 2;
    }
    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const long long n_genes = std::atoll(argv[3]);
    const long long n_cells = std::atoll(argv[4]);
    const std::string partial_path = output_path + ".partial";

    try {
        long long output_entries = aggregate_cells(
            input_path, n_genes, [&](long long, long long column, long long) {
                if (column < 1 || column > n_cells) {
                    throw std::runtime_error(
                        "MatrixMarket column is outside the declared cell range");
                }
            });

        BufferedGzipWriter writer(partial_path);
        writer.text("%%MatrixMarket matrix coordinate integer general\n");
        writer.text("% Aggregated from USA counts using gene-level S+A.\n");
        writer.text(std::to_string(n_genes) + " " + std::to_string(n_cells) + " " +
                    std::to_string(output_entries) + "\n");
        aggregate_cells(input_path, n_genes, [&](long long row, long long column, long long value) {
            writer.entry(row, column, value);
        });
        writer.flush();
        std::filesystem::rename(partial_path, output_path);
        std::cerr << "Output entries: " << output_entries << "\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
