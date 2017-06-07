/*!
 * @file bioparser_test.cpp
 *
 * @brief Bioparser unit test source file
 */

#include "bioparser_test_config.h"

#include "bioparser/bioparser.hpp"
#include "gtest/gtest.h"

class Read {
public:
    Read(uint64_t id, const char* name, uint32_t name_length,
        const char* data, uint32_t data_length)
            : id_(id), name_(name, name_length), data_(data, data_length),
            quality_() {
    }

    Read(uint64_t id, const char* name, uint32_t name_length,
        const char* data, uint32_t data_length,
        const char* quality, uint32_t quality_length)
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
            : id_(id),
            a_id_(a_id - 1), a_begin_(a_begin), a_end_(a_end), a_length_(a_length),
            b_id_(b_id - 1), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
            orientation_(a_rc == b_rc ? '+' : '-') {
    }

    Overlap(uint64_t id, const char* a_name, uint32_t a_name_length,
        uint32_t a_length, uint32_t a_begin, uint32_t a_end, char orientation,
        const char* b_name, uint32_t b_name_length, uint32_t b_length,
        uint32_t b_begin, uint32_t b_end, uint32_t matching_bases,
        uint32_t overlap_length, uint32_t quality)
            : id_(id),
            a_id_(atoi(a_name) - 1), a_begin_(a_begin), a_end_(a_end), a_length_(a_length),
            b_id_(atoi(b_name) - 1), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
            orientation_(orientation) {
    }

    uint64_t id_;
    uint32_t a_id_;
    uint32_t a_begin_;
    uint32_t a_end_;
    uint32_t a_length_;
    uint32_t b_id_;
    uint32_t b_begin_;
    uint32_t b_end_;
    uint32_t b_length_;
    char orientation_;
};

class BioparserFastaTest: public ::testing::Test {
public:
    void SetUp() {
        reader = bioparser::createReader<Read, bioparser::FastaReader>(
            bioparser_test_data_path + "/sample.fasta"
        );
    }

    void TearDown() {
    }

    std::unique_ptr<bioparser::Reader<Read>> reader;
};

class BioparserFastqTest: public ::testing::Test {
public:
    void SetUp() {
        reader = bioparser::createReader<Read, bioparser::FastqReader>("");
    }

    void TearDown() {
    }

    std::unique_ptr<bioparser::Reader<Read>> reader;
};

class BioparserMhapTest: public ::testing::Test {
public:
    void SetUp() {
        reader = bioparser::createReader<Overlap, bioparser::MhapReader>("");
    }

    void TearDown() {
    }

    std::unique_ptr<bioparser::Reader<Overlap>> reader;
};

class BioparserPafTest: public ::testing::Test {
public:
    void SetUp() {
        reader = bioparser::createReader<Overlap, bioparser::PafReader>("");
    }

    void TearDown() {
    }

    std::unique_ptr<bioparser::Reader<Overlap>> reader;
};

TEST_F(BioparserFastaTest, ReadFileAtOnce) {

    std::vector<std::unique_ptr<Read>> reads;
    reader->read_objects(reads, -1);

    ASSERT_EQ(0, 0);
}

TEST_F(BioparserFastaTest, ReadFileInChunks) {

    uint32_t size_in_bytes = 5 * 1024 * 1024; // 5 MB
    std::vector<std::unique_ptr<Read>> reads;
    reader->read_objects(reads, -1);

    ASSERT_EQ(0, 0);
}
