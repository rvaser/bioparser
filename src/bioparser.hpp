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

namespace bioparser {

constexpr uint32_t kSmallBufferSize = 4 * 1024; // 4 kB
constexpr uint32_t kMediumBufferSize = 5 * 1024 * 1024; // 5 MB
constexpr uint32_t kLargeBufferSize = 500 * 1024 * 1024; // 500 MB

/*!
 * @brief Reader class
 */
template<class T>
class Reader;

template<class T, template<class T_> class U>
std::unique_ptr<Reader<T>> createReader(const std::string& path);

/*!
 * @brief Reader class specializations
 */
template<class T>
class FastaReader; // inherits Reader

template<class T>
class FastqReader; // inherits Reader

template<class T>
class MhapReader; // inherits Reader

template<class T>
class PafReader; // inherits Reader

// Implementation of all defined classes and methods above
template<class T>
class Reader {
    public:
        virtual ~Reader() {}
        virtual bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) = 0;
        bool read_objects(std::vector<std::shared_ptr<T>>& dst, uint64_t max_bytes);

    protected:
        Reader(FILE* input_file)
                : input_file_(input_file, fclose), buffer_(kSmallBufferSize, 0),
                num_objects_read_(0) {
        }
        Reader(const Reader&) = delete;
        const Reader& operator=(const Reader&) = delete;

        std::unique_ptr<FILE, int(*)(FILE*)> input_file_;
        std::vector<char> buffer_;
        uint64_t num_objects_read_;
};

template<class T>
bool Reader<T>::read_objects(std::vector<std::shared_ptr<T>>& dst, uint64_t max_bytes) {

    std::vector<std::unique_ptr<T>> tmp;
    auto ret = read_objects(tmp, max_bytes);

    for (auto& it: tmp) {
        dst.push_back(std::move(it));
    }
    return ret;
}

template<class T, template<class T_> class U>
std::unique_ptr<Reader<T>> createReader(const std::string& path) {

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file != nullptr && "Unable to open file");

    return std::unique_ptr<Reader<T>>(new U<T>(input_file));
}

template<class T>
class FastaReader: public Reader<T> {
    public:
        ~FastaReader() {}
        bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
        friend std::unique_ptr<Reader<T>> createReader<T, FastaReader>(const std::string& path);

    private:
        FastaReader(FILE* input_file)
                : Reader<T>(input_file), large_buffer_(kMediumBufferSize, 0) {
        }
        FastaReader(const FastaReader&) = delete;
        const FastaReader& operator=(const FastaReader&) = delete;

        std::vector<char> large_buffer_;
};

template<class T>
bool FastaReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string name(kSmallBufferSize, 0);
    uint32_t name_length = 0;
    bool is_name = true;

    char* data = &this->large_buffer_[0];
    uint32_t data_length = 0;

    // unique_ptr<FILE> to FILE*
    auto input_file = this->input_file_.get();
    bool is_end = feof(input_file);

    while (!is_end) {

        uint64_t read_bytes = fread(this->buffer_.data(), sizeof(char),
            kSmallBufferSize, input_file);
        is_end = feof(input_file);

        total_bytes += read_bytes;
        if (max_bytes != 0 && total_bytes > max_bytes) {
            fseek(input_file, -(current_bytes + read_bytes), SEEK_CUR);
            status = true;
            break;
        }

        for (uint32_t i = 0; i < read_bytes; ++i) {
            auto c = this->buffer_[i];

            if (!is_name && (c == '>' || (is_end && i == read_bytes - 1))) {

                while (name_length > 0 && isspace(name[name_length - 1])) {
                    --name_length;
                }

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    name.c_str(), name_length, (const char*) data, data_length)));

                this->num_objects_read_ += 1;
                current_bytes = 0;

                name_length = 0;
                is_name = true;
                data_length = 0;
            }

            if (is_name) {
                if (c == '\n') {
                    is_name = false;
                } else if (name_length == kSmallBufferSize) {
                    continue;
                } else if (!(name_length == 0 && (c == '>' || isspace(c)))) {
                    name[name_length++] = c;
                }
            } else {
                if (!isspace(c)) {
                    data[data_length++] = c;
                    if (data_length >= this->large_buffer_.size()) {
                        this->large_buffer_.resize(kLargeBufferSize, 0);
                        data = &this->large_buffer_[0];
                    }
                }
            }

            ++current_bytes;
        }
    }

    return status;
}

template<class T>
class FastqReader: public Reader<T> {
    public:
        ~FastqReader() {}
        bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
        friend std::unique_ptr<Reader<T>> createReader<T, FastqReader>(const std::string& path);

    private:
        FastqReader(FILE* input_file)
                : Reader<T>(input_file), large_buffer_1_(kMediumBufferSize, 0),
                large_buffer_2_(kMediumBufferSize, 0) {
        }
        FastqReader(const FastqReader&) = delete;
        const FastqReader& operator=(const FastqReader&) = delete;

        std::vector<char> large_buffer_1_;
        std::vector<char> large_buffer_2_;
};

template<class T>
bool FastqReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string name(kSmallBufferSize, 0);
    uint32_t name_length = 0;

    char* data = &this->large_buffer_1_[0];
    uint32_t data_length = 0;

    char* quality = &this->large_buffer_2_[0];
    uint32_t quality_length = 0;

    // unique_ptr<FILE> to FILE*
    auto input_file = this->input_file_.get();
    bool is_end = feof(input_file);

    bool is_valid = false;
    uint32_t line_number = 0;

    while (!is_end) {

        uint64_t read_bytes = fread(this->buffer_.data(), sizeof(char),
            kSmallBufferSize, input_file);
        is_end = feof(input_file);

        total_bytes += read_bytes;
        if (max_bytes != 0 && total_bytes > max_bytes) {
            fseek(input_file, -(current_bytes + read_bytes), SEEK_CUR);
            status = true;
            break;
        }

        for (uint32_t i = 0; i < read_bytes; ++i) {
            auto c = this->buffer_[i];

            if (c == '\n') {
                line_number = (line_number + 1) % 4;
                if (line_number == 0) {
                    is_valid = true;
                }
            } else {
                switch (line_number) {
                    case 0:
                        if (name_length < kSmallBufferSize) {
                            if (!(name_length == 0 && (c == '@' || isspace(c)))) {
                                name[name_length++] = c;
                            }
                        }
                        break;
                    case 1:
                        data[data_length++] = c;
                        if (data_length >= this->large_buffer_1_.size()) {
                            this->large_buffer_1_.resize(kLargeBufferSize, 0);
                            data = &this->large_buffer_1_[0];
                            this->large_buffer_2_.resize(kLargeBufferSize, 0);
                            quality = &this->large_buffer_2_[0];
                        }
                        break;
                    case 2:
                        // comment line starting with '+'
                        // do nothing
                        break;
                    case 3:
                        quality[quality_length++] = c;
                        break;
                    default:
                        // never reaches this case
                        break;
                }
            }

            ++current_bytes;

            if (is_valid) {

                while (name_length > 0 && isspace(name[name_length - 1])) {
                    --name_length;
                }
                while (data_length > 0 && isspace(data[data_length - 1])) {
                    --data_length;
                }
                while (quality_length > 0 && isspace(quality[quality_length - 1])) {
                    --quality_length;
                }

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    name.c_str(), name_length, (const char*) data, data_length,
                    (const char*) quality, quality_length)));

                this->num_objects_read_ += 1;
                current_bytes = 0;

                name_length = 0;
                data_length = 0;
                quality_length = 0;
                is_valid = false;
            }
        }
    }

    return status;
}

template<class T>
class MhapReader: public Reader<T> {
    public:
        ~MhapReader() {}
        bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
        friend std::unique_ptr<Reader<T>> createReader<T, MhapReader>(const std::string& path);

    private:
        MhapReader(FILE* input_file)
                : Reader<T>(input_file) {
        }
        MhapReader(const MhapReader&) = delete;
        const MhapReader& operator=(const MhapReader&) = delete;
};

template<class T>
bool MhapReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string line(kSmallBufferSize, 0);
    uint32_t line_length = 0;

    const uint32_t kMhapObjectLength = 12;
    uint32_t values_length = 0;

    // unique_ptr<FILE> to FILE*
    auto input_file = this->input_file_.get();
    bool is_end = feof(input_file);

    uint32_t a_id = 0, b_id = 0, minmers = 0, a_rc = 0, a_begin = 0, a_end = 0, a_length = 0,
        b_rc = 0, b_begin = 0, b_end = 0, b_length = 0;
    double error = 0;

    while (!is_end) {

        uint64_t read_bytes = fread(this->buffer_.data(), sizeof(char),
            kSmallBufferSize, input_file);
        is_end = feof(input_file);

        total_bytes += read_bytes;
        if (max_bytes != 0 && total_bytes > max_bytes) {
            fseek(input_file, -(current_bytes + read_bytes), SEEK_CUR);
            status = true;
            break;
        }

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {

            auto c = this->buffer_[i];
            ++current_bytes;

            if (c == '\n') {

                line[line_length] = 0;
                while (line_length > 0 && isspace(line[line_length - 1])) {
                    line[line_length - 1] = 0;
                    --line_length;
                }

                size_t start = 0, end = 0;
                while (true) {
                    end = line.find(values_length == kMhapObjectLength - 1 ? '\0' : ' ', start);
                    if (end == std::string::npos) {
                        break;
                    }
                    line[end] = 0;

                    switch (values_length) {
                        case 0:
                            a_id = atoi(&line[start]);
                            break;
                        case 1:
                            b_id = atoi(&line[start]);
                            break;
                        case 2:
                            error = atof(&line[start]);
                            break;
                        case 3:
                            minmers = atoi(&line[start]);
                            break;
                        case 4:
                            a_rc = atoi(&line[start]);
                            break;
                        case 5:
                            a_begin = atoi(&line[start]);
                            break;
                        case 6:
                            a_end = atoi(&line[start]);
                            break;
                        case 7:
                            a_length = atoi(&line[start]);
                            break;
                        case 8:
                            b_rc = atoi(&line[start]);
                            break;
                        case 9:
                            b_begin = atoi(&line[start]);
                            break;
                        case 10:
                            b_end = atoi(&line[start]);
                            break;
                        case 11:
                        default:
                            b_length = atoi(&line[start]);
                            break;
                    }
                    values_length++;
                    start = end + 1;
                }
                line_length = 0;
                assert(values_length == kMhapObjectLength && "Invalid format");

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    a_id, b_id, error, minmers, a_rc, a_begin, a_end, a_length,
                    b_rc, b_begin, b_end, b_length)));

                this->num_objects_read_ += 1;
                current_bytes = 0;

                values_length = 0;
            } else {
                line[line_length++] = c;
            }
        }
    }

    return status;
}

template<class T>
class PafReader: public Reader<T> {
    public:
        ~PafReader() {}
        bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
        friend std::unique_ptr<Reader<T>> createReader<T, PafReader>(const std::string& path);

    private:
        PafReader(FILE* input_file)
                : Reader<T>(input_file) {
        }
        PafReader(const PafReader&) = delete;
        const PafReader& operator=(const PafReader&) = delete;
};

template<class T>
bool PafReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string line(kSmallBufferSize, 0);
    uint32_t line_length = 0;

    const uint32_t kPafObjectLength = 12;
    uint32_t values_length = 0;

    // unique_ptr<FILE> to FILE*
    auto input_file = this->input_file_.get();
    bool is_end = feof(input_file);

    const char* a_name = nullptr;
    uint32_t a_name_length = 0;
    const char* b_name = nullptr;
    uint32_t b_name_length = 0;

    uint32_t a_length = 0, a_begin = 0, a_end = 0, b_length = 0, b_begin = 0,
        b_end = 0, matching_bases = 0, overlap_length = 0, quality = 0;
    char orientation = '\0';

    while (!is_end) {

        uint64_t read_bytes = fread(this->buffer_.data(), sizeof(char),
            kSmallBufferSize, input_file);
        is_end = feof(input_file);

        total_bytes += read_bytes;
        if (max_bytes != 0 && total_bytes > max_bytes) {
            fseek(input_file, -(current_bytes + read_bytes), SEEK_CUR);
            status = true;
            break;
        }

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {

            auto c = this->buffer_[i];
            ++current_bytes;

            if (c == '\n') {

                line[line_length] = 0;
                while (line_length > 0 && isspace(line[line_length - 1])) {
                    line[line_length - 1] = 0;
                    --line_length;
                }

                size_t start = 0, end = 0;
                while (true) {
                    end = line.find(values_length == kPafObjectLength - 1 ? '\0' : '\t', start);
                    if (end == std::string::npos) {
                        break;
                    }
                    line[end] = 0;

                    switch (values_length) {
                        case 0:
                            a_name = &line[start];
                            a_name_length = end - start;
                            break;
                        case 1:
                            a_length = atoi(&line[start]);
                            break;
                        case 2:
                            a_begin = atoi(&line[start]);
                            break;
                        case 3:
                            a_end = atoi(&line[start]);
                            break;
                        case 4:
                            orientation = line[start];
                            break;
                        case 5:
                            b_name = &line[start];
                            b_name_length = end - start;
                            break;
                        case 6:
                            b_length = atoi(&line[start]);
                            break;
                        case 7:
                            b_begin = atoi(&line[start]);
                            break;
                        case 8:
                            b_end = atoi(&line[start]);
                            break;
                        case 9:
                            matching_bases = atoi(&line[start]);
                            break;
                        case 10:
                            overlap_length = atoi(&line[start]);
                            break;
                        case 11:
                        default:
                            quality = atoi(&line[start]);
                            break;
                    }
                    values_length++;
                    if (values_length == kPafObjectLength) {
                        break;
                    }
                    start = end + 1;
                }
                line_length = 0;
                assert(values_length == kPafObjectLength && "Invalid format");

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    a_name, a_name_length, a_length, a_begin, a_end, orientation,
                    b_name, b_name_length, b_length, b_begin, b_end,
                    matching_bases, overlap_length, quality)));

                this->num_objects_read_ += 1;
                current_bytes = 0;

                values_length = 0;
                a_name_length = 0;
                b_name_length = 0;
            } else {
                line[line_length++] = c;
            }
        }
    }

    return status;
}

}
