#ifdef BIOPARSER_MAIN_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "bioparser.hpp"

using namespace BIOPARSER;

class Read {
public:
    Read(uint32_t id, const char* name, uint32_t name_length, const char* data,
        uint32_t data_length)
            : id_(id), name_(name, name_length), data_(data, data_length),
            quality_() {
    }

    Read(uint32_t id, const char* name, uint32_t name_length, const char* data,
        uint32_t data_length, const char* quality, uint32_t quality_length)
            : id_(id), name_(name, name_length), data_(data, data_length),
            quality_(quality, quality_length) {
    }

    uint32_t id_;
    std::string name_;
    std::string data_;
    std::string quality_;
};

class Overlap {
public:
    Overlap(uint32_t id, uint32_t a_id, uint32_t b_id, float error, uint32_t minmers,
        uint32_t a_rc, uint32_t a_start, uint32_t a_end, uint32_t a_length,
        uint32_t b_rc, uint32_t b_start, uint32_t b_end, uint32_t b_length)
            : id_(id), a_id_(a_id), b_id_(b_id), error_(error), minmers_(minmers),
            a_rc_(a_rc), a_start_(a_start), a_end_(a_end), a_length_(a_length),
            b_rc_(b_rc), b_start_(b_start), b_end_(b_end), b_length_(b_length) {
    }
    Overlap(uint32_t id, const std::vector<double>& values)
            : id_(id), a_id_(values[0]), b_id_(values[1]), error_(values[2]), minmers_(values[3]),
            a_rc_(values[4]), a_start_(values[5]), a_end_(values[6]), a_length_(values[7]),
            b_rc_(values[8]), b_start_(values[9]), b_end_(values[10]), b_length_(values[11]) {
    }

    uint32_t id_;
    uint32_t a_id_;
    uint32_t b_id_;
    float error_;
    uint32_t minmers_;
    uint32_t a_rc_;
    uint32_t a_start_;
    uint32_t a_end_;
    uint32_t a_length_;
    uint32_t b_rc_;
    uint32_t b_start_;
    uint32_t b_end_;
    uint32_t b_length_;
};

int main(int argc, char** argv) {

    timeval start, stop;

    gettimeofday(&start, nullptr);

    auto reader = createFastaReader<Read>(argv[1]);
    std::vector<std::unique_ptr<Read>> test;
    reader->read_objects(test, 10000000000000);

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "FASTA: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    gettimeofday(&start, nullptr);

    auto qreader = createFastqReader<Read>(argv[2]);
    std::vector<std::unique_ptr<Read>> test2;
    qreader->read_objects(test2, 10000000000000);

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "FASTQ: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    gettimeofday(&start, nullptr);

    auto mreader = createMhapReader<Overlap>(argv[3]);
    std::vector<std::unique_ptr<Overlap>> test3;
    mreader->read_objects(test3, 10000000000000);

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "MHAP: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    return 0;
}

#endif
