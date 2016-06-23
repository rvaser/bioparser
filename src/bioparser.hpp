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

constexpr uint32_t kSmallBufferSize = 4 * 1024; // 4 kB
constexpr uint32_t kLargeBufferSize = 500 * 1024 * 1024; // 500 MB

/*!
 * @brief Reader class
 */
template<class T>
class Reader {
public:
    ~Reader() {}
    virtual bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) = 0;

protected:
    Reader(FILE* input_file) :
            input_file_(input_file, fclose), buffer_(kSmallBufferSize, 0),
            num_objects_read_(0) {
    }
    Reader(const Reader&) = delete;
    const Reader& operator=(const Reader&) = delete;

    std::unique_ptr<FILE, int(*)(FILE*)> input_file_;
    std::vector<char> buffer_;
    uint32_t num_objects_read_;
};

/*!
 * @brief Reader class specializations
 */
template<class T>
class FastaReader; // inherits Reader

template<class T>
std::unique_ptr<Reader<T>> createFastaReader(const std::string& path);

template<class T>
class FastqReader; // inherits Reader

template<class T>
std::unique_ptr<Reader<T>> createFastqReader(const std::string& path);

template<class T>
class MhapReader; // inherits Reader

template<class T>
std::unique_ptr<Reader<T>> createMhapReader(const std::string& path);

// taken (and modified a bit) from http://www.leapsecond.com/tools/fast_atof.c
static double fast_atof(const char* p);

/*!
 * @brief Writer class
 */
class Writer {
};

// Implementation of all defined (Reader/Writer) methods above
template<class T>
class FastaReader: public Reader<T> {
public:
    ~FastaReader() {}
    bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
    friend std::unique_ptr<Reader<T>> createFastaReader<T>(const std::string& path);

protected:
    FastaReader(FILE* input_file) :
        Reader<T>(input_file) {
    }
    FastaReader(const FastaReader&) = delete;
    const FastaReader& operator=(const FastaReader&) = delete;
};

template<class T>
bool FastaReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string name(kSmallBufferSize, 0);
    uint32_t name_length = 0;
    bool is_name = true;

    std::string data(kLargeBufferSize, 0);
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

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {
            auto c = this->buffer_[i];

            if (!is_name && (c == '>' || (is_end && i == read_bytes - 1))) {

                while (name_length > 0 && isspace(name[name_length - 1])) {
                    --name_length;
                }
                while (data_length > 0 && isspace(data[data_length - 1])) {
                    --data_length;
                }

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    name.c_str(), name_length, data.c_str(), data_length)));

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
                data[data_length++] = c;
            }
        }

        current_bytes += i;
    }

    return status;
}

template<class T>
std::unique_ptr<Reader<T>> createFastaReader(const std::string& path) {

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file != nullptr && "Unable to open file");

    return std::unique_ptr<Reader<T>>(new FastaReader<T>(input_file));
}

template<class T>
class FastqReader: public Reader<T> {
public:
    ~FastqReader() {}
    bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
    friend std::unique_ptr<Reader<T>> createFastqReader<T>(const std::string& path);

protected:
    FastqReader(FILE* input_file) :
        Reader<T>(input_file) {
    }
    FastqReader(const FastqReader&) = delete;
    const FastqReader& operator=(const FastqReader&) = delete;
};

template<class T>
bool FastqReader<T>::read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes) {

    bool status = false;
    uint64_t current_bytes = 0;
    uint64_t total_bytes = 0;

    std::string name(kSmallBufferSize, 0);
    uint32_t name_length = 0;

    std::string data(kLargeBufferSize, 0);
    uint32_t data_length = 0;

    std::string quality(kLargeBufferSize, 0);
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

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {
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
                    name.c_str(), name_length, data.c_str(), data_length,
                    quality.c_str(), quality_length)));

                this->num_objects_read_ += 1;

                current_bytes = 0;
                name_length = 0;
                data_length = 0;
                quality_length = 0;
                is_valid = false;
            }
        }

        current_bytes += i;
    }

    return status;
}

template<class T>
std::unique_ptr<Reader<T>> createFastqReader(const std::string& path) {

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file != nullptr && "Unable to open file");

    return std::unique_ptr<Reader<T>>(new FastqReader<T>(input_file));
}

template<class T>
class MhapReader: public Reader<T> {
public:
    ~MhapReader() {}
    bool read_objects(std::vector<std::unique_ptr<T>>& dst, uint64_t max_bytes);
    friend std::unique_ptr<Reader<T>> createMhapReader<T>(const std::string& path);

protected:
    MhapReader(FILE* input_file) :
        Reader<T>(input_file) {
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
    std::vector<double> values(kMhapObjectLength, 0);
    uint32_t values_length = 0;

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

        uint32_t i = 0;
        for (; i < read_bytes; ++i) {

            auto c = this->buffer_[i];

            if (c == '\n') {

                line[line_length] = 0;
                while (line_length > 0 && isspace(line[line_length - 1])) {
                    line[line_length - 1] = 0;
                    --line_length;
                }

                int32_t start = -1, end = 0;
                while (true) {
                    end = line.find(values_length == kMhapObjectLength - 1 ? '\0' : ' ', start + 1);
                    if (end == std::string::npos) {
                        break;
                    }
                    line[end] = 0;
                    values[values_length++] = fast_atof(&line[start + 1]);
                    start = end;
                }
                line_length = 0;

                assert(values_length == kMhapObjectLength && "Invalid format");

                dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    (const double*) values.data(), values_length)));

                /*dst.emplace_back(std::unique_ptr<T>(new T(this->num_objects_read_,
                    (uint32_t) values[0], (uint32_t) values[1], // A id, B id
                               values[2], (uint32_t) values[3], // % error, # shared min-mers
                    (uint32_t) values[4], (uint32_t) values[5], // A fwd/rc, A start
                    (uint32_t) values[6], (uint32_t) values[7], // A end, A length
                    (uint32_t) values[8], (uint32_t) values[9], // B fwd/rc, B start
                    (uint32_t) values[10], (uint32_t) values[11]))); // B end, B length
*/
                this->num_objects_read_ += 1;
                values_length = 0;
            } else {
                line[line_length++] = c;
            }
        }

        current_bytes += i;
    }

    return status;
}

template<class T>
std::unique_ptr<Reader<T>> createMhapReader(const std::string& path) {

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file != nullptr && "Unable to open file");

    return std::unique_ptr<Reader<T>>(new MhapReader<T>(input_file));
}

static double fast_atof(const char* p) {
    int frac;
    double sign, value, scale;

    auto valid_digit = [](char c) -> bool {
        return (c >= '0' && c <= '9');
    };

    // Skip leading white space, if any.
    while (isspace(*p)) {
        ++p;
    }

    // Get sign, if any.
    sign = 1.0;
    if (*p == '-') {
        sign = -1.0;
        ++p;
    } else if (*p == '+') {
        ++p;
    }

    // Get digits before decimal point or exponent, if any.
    for (value = 0.0; valid_digit(*p); ++p) {
        value = value * 10.0 + (*p - '0');
    }

    // Get digits after decimal point, if any.
    if (*p == '.') {
        double pow10 = 10.0;
        ++p;
        while (valid_digit(*p)) {
            value += (*p - '0') / pow10;
            pow10 *= 10.0;
            ++p;
        }
    }

    // Handle exponent, if any.
    frac = 0;
    scale = 1.0;
    if ((*p == 'e') || (*p == 'E')) {
        unsigned int expon;

        // Get sign of exponent, if any.
        ++p;
        if (*p == '-') {
            frac = 1;
            ++p;
        } else if (*p == '+') {
            ++p;
        }

        // Get digits of exponent, if any.
        for (expon = 0; valid_digit(*p); ++p) {
            expon = expon * 10 + (*p - '0');
        }
        if (expon > 308) expon = 308;

        // Calculate scaling factor.
        while (expon >= 50) { scale *= 1E50; expon -= 50; }
        while (expon >=  8) { scale *= 1E8;  expon -=  8; }
        while (expon >   0) { scale *= 10.0; expon -=  1; }
    }

    // Return signed and scaled floating point result.
    return sign * (frac ? (value / scale) : (value * scale));
}

}
