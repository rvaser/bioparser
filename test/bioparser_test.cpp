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
    Read(uint64_t id_, const char* name_, uint32_t name_length_,
        const char* data_, uint32_t data_length_)
            : id(id_), name(name_, name_length_), data(data_, data_length_),
            quality() {
    }

    Read(uint64_t id_, const char* name_, uint32_t name_length_,
        const char* data_, uint32_t data_length_,
        const char* quality_, uint32_t quality_length_)
            : id(id_), name(name_, name_length_), data(data_, data_length_),
            quality(quality_, quality_length_) {
    }

    ~Read() {}

    uint64_t id;
    std::string name;
    std::string data;
    std::string quality;
};

void reads_summary(uint32_t& name_size, uint32_t& data_size, uint32_t& quality_size,
    const std::vector<std::unique_ptr<Read>>& reads) {

    name_size = 0;
    data_size = 0;
    quality_size = 0;
    for (const auto& it: reads) {
        name_size += it->name.size();
        data_size += it->data.size();
        quality_size += it->quality.size();
    }
}

class Overlap {
public:
    Overlap(uint64_t id_, uint32_t a_id_, uint32_t b_id_, double error, uint32_t minmers,
        uint32_t a_rc, uint32_t a_begin_, uint32_t a_end_, uint32_t a_length_,
        uint32_t b_rc, uint32_t b_begin_, uint32_t b_end_, uint32_t b_length_)
            : id(id_),
            a_id(a_id_ - 1), a_begin(a_begin_), a_end(a_end_), a_length(a_length_),
            b_id(b_id_ - 1), b_begin(b_begin_), b_end(b_end_), b_length(b_length_),
            orientation(a_rc == b_rc ? '+' : '-') {
    }

    Overlap(uint64_t id_, const char* a_name, uint32_t a_name_length,
        uint32_t a_length_, uint32_t a_begin_, uint32_t a_end_, char orientation_,
        const char* b_name, uint32_t b_name_length, uint32_t b_length_,
        uint32_t b_begin_, uint32_t b_end_, uint32_t matching_bases,
        uint32_t overlap_length, uint32_t quality)
            : id(id_),
            a_id(atoi(a_name) - 1), a_begin(a_begin_), a_end(a_end_), a_length(a_length_),
            b_id(atoi(b_name) - 1), b_begin(b_begin_), b_end(b_end_), b_length(b_length_),
            orientation(orientation_) {
    }

    ~Overlap() {}

    uint64_t id;
    uint32_t a_id;
    uint32_t a_begin;
    uint32_t a_end;
    uint32_t a_length;
    uint32_t b_id;
    uint32_t b_begin;
    uint32_t b_end;
    uint32_t b_length;
    char orientation;
};

class BioparserFastaTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        reader = bioparser::createReader<Read, bioparser::FastaReader>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Reader<Read>> reader;
};

class BioparserFastqTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        reader = bioparser::createReader<Read, bioparser::FastqReader>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Reader<Read>> reader;
};

class BioparserMhapTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        reader = bioparser::createReader<Overlap, bioparser::MhapReader>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Reader<Overlap>> reader;
};

class BioparserPafTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        reader = bioparser::createReader<Overlap, bioparser::PafReader>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Reader<Overlap>> reader;
};

TEST(BioparserTest, CreateReaderError) {
    EXPECT_DEATH((bioparser::createReader<Read, bioparser::FastaReader>("")),
        "bioparser::createReader error: unable to open file");
}

TEST_F(BioparserFastaTest, ReadWhole) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    std::vector<std::unique_ptr<Read>> reads;
    reader->read_objects(reads, -1);

    uint32_t reads_name_size = 0, reads_data_size = 0, reads_quality_size = 0;
    reads_summary(reads_name_size, reads_data_size, reads_quality_size, reads);

    EXPECT_EQ(14U, reads.size());
    EXPECT_EQ(75U, reads_name_size);
    EXPECT_EQ(109117U, reads_data_size);
    EXPECT_EQ(0U, reads_quality_size);
}

TEST_F(BioparserFastaTest, ReadInChunks) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    uint32_t size_in_bytes = 25 * 1024; // 25 kB
    std::vector<std::unique_ptr<Read>> reads;
    while (reader->read_objects(reads, size_in_bytes)) {
    }

    uint32_t reads_name_size = 0, reads_data_size = 0, reads_quality_size = 0;
    reads_summary(reads_name_size, reads_data_size, reads_quality_size, reads);

    EXPECT_EQ(14U, reads.size());
    EXPECT_EQ(75U, reads_name_size);
    EXPECT_EQ(109117U, reads_data_size);
    EXPECT_EQ(0U, reads_quality_size);
}

TEST_F(BioparserFastaTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.fastq");
    std::vector<std::unique_ptr<Read>> reads;

    EXPECT_DEATH(reader->read_objects(reads, -1),
        "bioparser::FastaReader error: invalid file format!");
}

TEST_F(BioparserFastaTest, ChunkSizeError) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    uint32_t size_in_bytes = 10 * 1024; // 10 kB
    std::vector<std::unique_ptr<Read>> reads;
    EXPECT_DEATH(reader->read_objects(reads, size_in_bytes),
        "bioparser::FastaReader error: too small chunk size!");
}

TEST_F(BioparserFastqTest, ReadWhole) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    std::vector<std::unique_ptr<Read>> reads;
    reader->read_objects(reads, -1);

    uint32_t reads_name_size = 0, reads_data_size = 0, reads_quality_size = 0;
    reads_summary(reads_name_size, reads_data_size, reads_quality_size, reads);

    EXPECT_EQ(13U, reads.size());
    EXPECT_EQ(17U, reads_name_size);
    EXPECT_EQ(108140U, reads_data_size);
    EXPECT_EQ(108140U, reads_quality_size);
}

TEST_F(BioparserFastqTest, ReadInChunks) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    uint32_t size_in_bytes = 50 * 1024; // 50 kB
    std::vector<std::unique_ptr<Read>> reads;
    while (reader->read_objects(reads, size_in_bytes)) {
    }

    uint32_t reads_name_size = 0, reads_data_size = 0, reads_quality_size = 0;
    reads_summary(reads_name_size, reads_data_size, reads_quality_size, reads);

    EXPECT_EQ(13U, reads.size());
    EXPECT_EQ(17U, reads_name_size);
    EXPECT_EQ(108140U, reads_data_size);
    EXPECT_EQ(108140U, reads_quality_size);
}

TEST_F(BioparserFastqTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    std::vector<std::unique_ptr<Read>> reads;

    EXPECT_DEATH(reader->read_objects(reads, -1),
        "bioparser::FastqReader error: invalid file format!");
}

TEST_F(BioparserFastqTest, ChunkSizeError) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    uint32_t size_in_bytes = 10 * 1024; // 10 kB
    std::vector<std::unique_ptr<Read>> reads;
    EXPECT_DEATH(reader->read_objects(reads, size_in_bytes),
        "bioparser::FastqReader error: too small chunk size!");
}

TEST_F(BioparserMhapTest, ReadWhole) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    std::vector<std::unique_ptr<Overlap>> overlaps;
    reader->read_objects(overlaps, -1);

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserMhapTest, ReadFileInChunks) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    uint32_t size_in_bytes = 4 * 1024; // 1 kB
    std::vector<std::unique_ptr<Overlap>> overlaps;
    while (reader->read_objects(overlaps, size_in_bytes)) {
    }

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserMhapTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.paf");

    std::vector<std::unique_ptr<Overlap>> overlaps;

    EXPECT_DEATH(reader->read_objects(overlaps, -1),
        "bioparser::MhapReader error: invalid file format!");
}

TEST_F(BioparserPafTest, ReadFileAtOnce) {

    SetUp(bioparser_test_data_path + "sample.paf");

    std::vector<std::unique_ptr<Overlap>> overlaps;
    reader->read_objects(overlaps, -1);

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserPafTest, ReadFileInChunks) {

    SetUp(bioparser_test_data_path + "sample.paf");

    uint32_t size_in_bytes = 4 * 1024; // 1 kB
    std::vector<std::unique_ptr<Overlap>> overlaps;
    while (reader->read_objects(overlaps, size_in_bytes)) {
    }

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserPafTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    std::vector<std::unique_ptr<Overlap>> overlaps;

    EXPECT_DEATH(reader->read_objects(overlaps, -1),
        "bioparser::PafReader error: invalid file format!");
}
