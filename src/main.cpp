#ifdef BIOPARSER_MAIN_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "bioparser.hpp"

using namespace bioparser;

class Read {
    public:
        Read(uint64_t id, const char* name, uint32_t name_length, const char* data,
            uint32_t data_length)
                : id_(id), name_(name, name_length), data_(data, data_length), quality_() {
        }

        Read(uint64_t id, const char* name, uint32_t name_length, const char* data,
            uint32_t data_length, const char* quality, uint32_t quality_length)
                : id_(id), name_(name, name_length), data_(data, data_length),
                quality_(quality, quality_length) {
        }

        uint64_t id_;
        std::string name_;
        std::string data_;
        std::string quality_;
};

class Overlap {
    public:
        Overlap(uint64_t id, uint32_t a_id, uint32_t b_id, double error, uint32_t minmers,
            uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
            uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
                : id_(id), a_id_(a_id), b_id_(b_id), error_(error), minmers_(minmers),
                a_rc_(a_rc), a_begin_(a_begin), a_end_(a_end), a_length_(a_length),
                b_rc_(b_rc), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
                quality_(), a_name_(), b_name_() {
        }

        Overlap(uint64_t id, const char* a_name, uint32_t a_name_length, uint32_t a_length, uint32_t a_begin, uint32_t a_end, char orientation,
            const char* b_name, uint32_t b_name_length, uint32_t b_length, uint32_t b_begin, uint32_t b_end,
            uint32_t matching_bases, uint32_t overlap_length, uint32_t quality)
                : id_(id), a_id_(atoi(a_name)), b_id_(atoi(b_name)), error_(matching_bases / (double) overlap_length), minmers_(),
                a_rc_(), a_begin_(a_begin), a_end_(a_end), a_length_(a_length),
                b_rc_(orientation == '-' ? 1 : 0), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
                quality_(quality), a_name_(a_name, a_name_length), b_name_(b_name, b_name_length) {
        }

        uint64_t id_;
        uint32_t a_id_;
        uint32_t b_id_;
        double error_;
        uint32_t minmers_;
        uint32_t a_rc_;
        uint32_t a_begin_;
        uint32_t a_end_;
        uint32_t a_length_;
        uint32_t b_rc_;
        uint32_t b_begin_;
        uint32_t b_end_;
        uint32_t b_length_;
        uint32_t quality_;
        std::string a_name_;
        std::string b_name_;
};

int main(int argc, char** argv) {

    timeval start, stop;

    gettimeofday(&start, nullptr);

    auto reader = createReader<Read, FastaReader>(argv[1]);
    std::vector<std::unique_ptr<Read>> test;

    uint32_t size_in_bytes = 5 * 1024 * 1024; // 5 MB
    //uint32_t size_in_bytes = -1; // everything
    while (true) {
        auto status = reader->read_objects(test, size_in_bytes);
        if (!status) {
            break;
        }
    }

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "FASTA: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    /*for (const auto& it: test) {
        fprintf(stdout, ">%s\n%s\n", it->name_.c_str(), it->data_.c_str());
    }*/

    gettimeofday(&start, nullptr);

    auto qreader = createReader<Read, FastqReader>(argv[2]);
    std::vector<std::unique_ptr<Read>> test2;
    while (true) {
        auto status = qreader->read_objects(test2, size_in_bytes);
        if (!status) {
            break;
        }
    }

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "FASTQ: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    /*for (const auto& it: test2) {
        fprintf(stdout, "@%s\n%s\n+\n%s\n", it->name_.c_str(), it->data_.c_str(), it->quality_.c_str());
    }*/

    gettimeofday(&start, nullptr);

    auto mreader = createReader<Overlap, MhapReader>(argv[3]);
    std::vector<std::unique_ptr<Overlap>> test3;
    while (true) {
        auto status = mreader->read_objects(test3, size_in_bytes);
        if (!status) {
            break;
        }
    }

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "MHAP: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    /*for (const auto& it: test3) {
        fprintf(stdout, "%d %d %g %d %d %d %d %d %d %d %d %d\n", it->a_id_, it->b_id_, it->error_, it->minmers_,
            it->a_rc_, it->a_begin_, it->a_end_, it->a_length_,
            it->b_rc_, it->b_begin_, it->b_end_, it->b_length_);
    }*/

    gettimeofday(&start, nullptr);

    auto preader = createReader<Overlap, PafReader>(argv[4]);
    std::vector<std::unique_ptr<Overlap>> test4;
    while (true) {
        auto status = preader->read_objects(test4, size_in_bytes);
        if (!status) {
            break;
        }
    }

    gettimeofday(&stop, nullptr);
    fprintf(stderr, "PAF: %.5lf\n", (((stop.tv_sec - start.tv_sec) * 1000000L + stop.tv_usec)
        - start.tv_usec) / (double) 1000000);

    /*for (const auto& it: test4) {
        uint32_t overlap_length = std::max(it->a_end_ - it->a_begin_, it->b_end_ - it->b_begin_);
        fprintf(stdout, "%s %d %d %d %c %s %d %d %d %d %d %d\n", it->a_name_.c_str(), it->a_length_, it->a_begin_, it->a_end_,
            (it->b_rc_ == 1 ? '-' : '+'), it->b_name_.c_str(), it->b_length_, it->b_begin_, it->b_end_,
            (uint32_t) (it->error_ * overlap_length + 0.499), overlap_length, it->quality_);
    }*/

    return 0;
}

#endif
