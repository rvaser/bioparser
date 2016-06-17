/*!
 * @file bioparser.hpp
 *
 * @brief Bioparser header file
 */

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <memory>
#include <string>
#include <vector>

namespace BIOPARSER {

constexpr uint32_t kBufferSize = 4 * 1024 * 1024; // 4kB

enum class FileType {
    kFASTA,
    kFASTQ,
    kMHAP,
};

/*!
 * @brief Reader class
 */
template<class T>
class Reader;

template<class T>
std::unique_ptr<Reader<T>> createReader(const std::string& path, FileType type);

template<class T>
class Reader {
public:
    ~Reader() {}

    bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);

    template<class T_>
    friend std::unique_ptr<Reader<T_>> createReader(const std::string& path,
        FileType type);

private:
    Reader(FILE* input_file, FileType type)
            : input_file_(input_file, fclose), buffer_(kBufferSize, '0'),
            num_objects_read_(0), type_(type) {
    }
    Reader(const Reader&) = delete;
    const Reader& operator=(const Reader&) = delete;

    bool read_FASTA_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
    bool read_FASTQ_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
    bool read_MHAP_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);

    std::unique_ptr<FILE, int(*)(FILE*)> input_file_;
    std::vector<char> buffer_;
    uint32_t num_objects_read_;
    FileType type_;
};

/*!
 * @brief Writer class
 */
class Writer {
};

// Implementation of all defined (Reader/Writer) methods above
template<class T>
std::unique_ptr<Reader<T>> createReader(const std::string& path, FileType type) {

    assert((type == FileType::kFASTA || type == FileType::kFASTQ ||
        type == FileType::kMHAP) && "Unsupported file format");

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file != nullptr && "Unable to open file");

    return std::unique_ptr<Reader<T>>(new Reader<T>(input_file, type));
}

template<class T>
bool Reader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {
    switch (type_) {
        case FileType::kFASTA:
            return read_FASTA_objects(dst, max_bytes);
        case FileType::kFASTQ:
            return read_FASTQ_objects(dst, max_bytes);
        case FileType::kMHAP:
            return read_MHAP_objects(dst, max_bytes);
        default:
            // never reaches this case due to asserts in createReader
            return false;
    }
}

template<class T>
bool Reader<T>::read_FASTA_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string name, data;
    name.reserve(kBufferSize);
    data.reserve(kBufferSize);
    bool is_name = true;

    // unique_ptr<FILE> to FILE*
    auto input_file = input_file_.get();
    bool is_end = feof(input_file);

    while (!is_end) {

        uint64_t read_bytes = fread(buffer_.data(), sizeof(char), kBufferSize, input_file);
        is_end = feof(input_file);

        total_bytes += read_bytes;
        if (max_bytes != 0 && total_bytes > max_bytes) {
            fseek(input_file, -(current_bytes + read_bytes), SEEK_CUR);
            status = true;
            break;
        }

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {
            auto c = buffer_[i];

            if (!is_name && (c == '>' || (is_end && i == read_bytes - 1))) {
                current_bytes = 0;
                is_name = true;
                dst.emplace_back(std::unique_ptr<T>(new T(name, data)));
                name.clear();
                data.clear();
            }

            if (is_name) {
                if (c == '\n') {
                    is_name = false;
                } else if (name.size() == kBufferSize) {
                    continue;
                } else if (!(name.size() == 0 && (c == '>' || isspace(c))) && c != '\r') {
                    name.push_back(c);
                }
            } else {
                data.push_back(c);
            }
        }

        current_bytes += i;
    }

    return status;
}

template<class T>
bool Reader<T>::read_FASTQ_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {
    return false;
}

template<class T>
bool Reader<T>::read_MHAP_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {
    return false;
}

}
