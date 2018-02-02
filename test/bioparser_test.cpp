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
    Read(const char* name, uint32_t name_length, const char* sequence,
        uint32_t sequence_length)
            : name_(name, name_length), sequence_(sequence, sequence_length),
            quality_() {
    }

    Read(const char* name, uint32_t name_length, const char* sequence,
        uint32_t sequence_length, const char* quality, uint32_t quality_length)
            : name_(name, name_length), sequence_(sequence, sequence_length),
            quality_(quality, quality_length) {
    }

    ~Read() {}

    std::string name_;
    std::string sequence_;
    std::string quality_;
};

void reads_summary(uint32_t& name_size, uint32_t& sequence_size, uint32_t& quality_size,
    const std::vector<std::unique_ptr<Read>>& reads) {

    name_size = 0;
    sequence_size = 0;
    quality_size = 0;
    for (const auto& it: reads) {
        name_size += it->name_.size();
        sequence_size += it->sequence_.size();
        quality_size += it->quality_.size();
    }
}

class Overlap {
public:
    Overlap(uint64_t a_id, uint64_t b_id, double error, uint32_t minmers,
        uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
        uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
            : a_name_(), a_id_(a_id - 1), a_begin_(a_begin), a_end_(a_end),
            a_length_(a_length), b_name_(), b_id_(b_id - 1), b_begin_(b_begin),
            b_end_(b_end), b_length_(b_length), orientation_(a_rc == b_rc ? '+' : '-'),
            error_(error), minmers_(minmers), matching_bases_(), overlap_length_(),
            mapping_quality_() {
    }

    Overlap(const char* a_name, uint32_t a_name_length, uint32_t a_length,
        uint32_t a_begin, uint32_t a_end, char orientation, const char* b_name,
        uint32_t b_name_length, uint32_t b_length, uint32_t b_begin,
        uint32_t b_end, uint32_t matching_bases, uint32_t overlap_length,
        uint32_t quality)
            : a_name_(a_name, a_name_length), a_id_(), a_begin_(a_begin),
            a_end_(a_end), a_length_(a_length), b_name_(b_name, b_name_length),
            b_id_(), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
            orientation_(orientation), error_(), minmers_(),
            matching_bases_(matching_bases), overlap_length_(overlap_length),
            mapping_quality_(quality) {
    }

    ~Overlap() {}

    std::string a_name_;
    uint64_t a_id_;
    uint32_t a_begin_;
    uint32_t a_end_;
    uint32_t a_length_;
    std::string b_name_;
    uint64_t b_id_;
    uint32_t b_begin_;
    uint32_t b_end_;
    uint32_t b_length_;
    char orientation_;
    uint32_t error_;
    uint32_t minmers_;
    uint32_t matching_bases_;
    uint32_t overlap_length_;
    uint32_t mapping_quality_;
};

class Alignment {
public:
    Alignment(const char* q_name, uint32_t q_name_length, uint32_t flag,
        const char* t_name, uint32_t t_name_length, uint32_t t_begin,
        uint32_t mapping_quality, const char* cigar, uint32_t cigar_length,
        const char* t_next_name, uint32_t t_next_name_length, uint32_t t_next_begin,
        uint32_t template_length, const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length)
            : q_name_(q_name, q_name_length), flag_(flag), t_name_(t_name, t_name_length),
            t_begin_(t_begin), mapping_quality_(mapping_quality),
            cigar_(cigar, cigar_length), t_next_name_(t_next_name, t_next_name_length),
            t_next_begin_(t_next_begin), template_length_(template_length),
            sequence_(sequence, sequence_length), quality_(quality, quality_length) {
    }

    ~Alignment() {}

    std::string q_name_;
    uint32_t flag_;
    std::string t_name_;
    uint32_t t_begin_;
    uint32_t mapping_quality_;
    std::string cigar_;
    std::string t_next_name_;
    uint32_t t_next_begin_;
    uint32_t template_length_;
    std::string sequence_;
    std::string quality_;
};

class BioparserFastaTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        parser = bioparser::createParser<bioparser::FastaParser, Read>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Parser<Read>> parser;
};

class BioparserFastqTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        parser = bioparser::createParser<bioparser::FastqParser, Read>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Parser<Read>> parser;
};

class BioparserMhapTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        parser = bioparser::createParser<bioparser::MhapParser, Overlap>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Parser<Overlap>> parser;
};

class BioparserPafTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        parser = bioparser::createParser<bioparser::PafParser, Overlap>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Parser<Overlap>> parser;
};

class BioparserSamTest: public ::testing::Test {
public:
    void SetUp(const std::string& file_name) {
        parser = bioparser::createParser<bioparser::SamParser, Alignment>(file_name);
    }

    void TearDown() {}

    std::unique_ptr<bioparser::Parser<Alignment>> parser;
};

TEST(BioparserTest, CreateParserError) {
    EXPECT_DEATH((bioparser::createParser<bioparser::FastaParser, Read>("")),
        ".bioparser::createParser. error: unable to open file !");
}

TEST_F(BioparserFastaTest, ParseWhole) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    std::vector<std::unique_ptr<Read>> reads;
    parser->parse_objects(reads, -1);

    uint32_t name_size = 0, sequence_size = 0, quality_size = 0;
    reads_summary(name_size, sequence_size, quality_size, reads);

    EXPECT_EQ(14U, reads.size());
    EXPECT_EQ(65U, name_size);
    EXPECT_EQ(109117U, sequence_size);
    EXPECT_EQ(0U, quality_size);
}

TEST_F(BioparserFastaTest, ParseInChunks) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    uint32_t size_in_bytes = 64 * 1024;
    std::vector<std::unique_ptr<Read>> reads;
    while (parser->parse_objects(reads, size_in_bytes)) {
    }

    uint32_t name_size = 0, sequence_size = 0, quality_size = 0;
    reads_summary(name_size, sequence_size, quality_size, reads);

    EXPECT_EQ(14U, reads.size());
    EXPECT_EQ(65U, name_size);
    EXPECT_EQ(109117U, sequence_size);
    EXPECT_EQ(0U, quality_size);
}

TEST_F(BioparserFastaTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.fastq");
    std::vector<std::unique_ptr<Read>> reads;

    EXPECT_DEATH(parser->parse_objects(reads, -1),
        ".bioparser::FastaParser. error: invalid file format!");
}

TEST_F(BioparserFastaTest, ChunkSizeError) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    uint32_t size_in_bytes = 10 * 1024;
    std::vector<std::unique_ptr<Read>> reads;
    EXPECT_DEATH(parser->parse_objects(reads, size_in_bytes),
        ".bioparser::FastaParser. error: too small chunk size!");
}

TEST_F(BioparserFastaTest, ParseAndReset) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    std::vector<std::unique_ptr<Read>> reads;
    parser->parse_objects(reads, -1);

    uint32_t num_reads = reads.size(), name_size = 0, sequence_size = 0,
        quality_size = 0;
    reads_summary(name_size, sequence_size, quality_size, reads);

    uint32_t size_in_bytes = 64 * 1024;
    reads.clear();
    parser->reset();
    while (parser->parse_objects(reads, size_in_bytes)) {
    }

    uint32_t num_reads_new = reads.size(), name_size_new = 0,
        sequence_size_new = 0, quality_size_new = 0;
    reads_summary(name_size_new, sequence_size_new, quality_size_new, reads);

    EXPECT_EQ(num_reads_new, num_reads);
    EXPECT_EQ(name_size_new, name_size);
    EXPECT_EQ(sequence_size_new, sequence_size);
    EXPECT_EQ(quality_size_new, quality_size);
}

TEST_F(BioparserFastqTest, ParseWhole) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    std::vector<std::unique_ptr<Read>> reads;
    parser->parse_objects(reads, -1);

    uint32_t name_size = 0, sequence_size = 0, quality_size = 0;
    reads_summary(name_size, sequence_size, quality_size, reads);

    EXPECT_EQ(13U, reads.size());
    EXPECT_EQ(17U, name_size);
    EXPECT_EQ(108140U, sequence_size);
    EXPECT_EQ(108140U, quality_size);
}

TEST_F(BioparserFastqTest, ParseInChunks) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    uint32_t size_in_bytes = 64 * 1024;
    std::vector<std::unique_ptr<Read>> reads;
    while (parser->parse_objects(reads, size_in_bytes)) {
    }

    uint32_t name_size = 0, sequence_size = 0, quality_size = 0;
    reads_summary(name_size, sequence_size, quality_size, reads);

    EXPECT_EQ(13U, reads.size());
    EXPECT_EQ(17U, name_size);
    EXPECT_EQ(108140U, sequence_size);
    EXPECT_EQ(108140U, quality_size);
}

TEST_F(BioparserFastqTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.fasta");

    std::vector<std::unique_ptr<Read>> reads;

    EXPECT_DEATH(parser->parse_objects(reads, -1),
        ".bioparser::FastqParser. error: invalid file format!");
}

TEST_F(BioparserFastqTest, ChunkSizeError) {

    SetUp(bioparser_test_data_path + "sample.fastq");

    uint32_t size_in_bytes = 10 * 1024;
    std::vector<std::unique_ptr<Read>> reads;
    EXPECT_DEATH(parser->parse_objects(reads, size_in_bytes),
        ".bioparser::FastqParser. error: too small chunk size!");
}

TEST_F(BioparserMhapTest, ParseWhole) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    std::vector<std::unique_ptr<Overlap>> overlaps;
    parser->parse_objects(overlaps, -1);

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserMhapTest, ParseInChunks) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    uint32_t size_in_bytes = 64 * 1024;
    std::vector<std::unique_ptr<Overlap>> overlaps;
    while (parser->parse_objects(overlaps, size_in_bytes)) {
    }

    EXPECT_EQ(150U, overlaps.size());
}

TEST_F(BioparserMhapTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.paf");

    std::vector<std::unique_ptr<Overlap>> overlaps;

    EXPECT_DEATH(parser->parse_objects(overlaps, -1),
        ".bioparser::MhapParser. error: invalid file format!");
}

TEST_F(BioparserPafTest, ParseWhole) {

    SetUp(bioparser_test_data_path + "sample.paf");

    std::vector<std::unique_ptr<Overlap>> overlaps;
    parser->parse_objects(overlaps, -1);

    EXPECT_EQ(500U, overlaps.size());
}

TEST_F(BioparserPafTest, ParseInChunks) {

    SetUp(bioparser_test_data_path + "sample.paf");

    uint32_t size_in_bytes = 64 * 1024;
    std::vector<std::unique_ptr<Overlap>> overlaps;
    while (parser->parse_objects(overlaps, size_in_bytes)) {
    }

    EXPECT_EQ(500U, overlaps.size());
}

TEST_F(BioparserPafTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.mhap");

    std::vector<std::unique_ptr<Overlap>> overlaps;

    EXPECT_DEATH(parser->parse_objects(overlaps, -1),
        ".bioparser::PafParser. error: invalid file format!");
}

TEST_F(BioparserSamTest, ParseWhole) {

    SetUp(bioparser_test_data_path + "sample.sam");

    std::vector<std::unique_ptr<Alignment>> alignments;
    parser->parse_objects(alignments, -1);

    EXPECT_EQ(48U, alignments.size());
}

TEST_F(BioparserSamTest, ParseInChunks) {

    SetUp(bioparser_test_data_path + "sample.sam");

    uint32_t size_in_bytes = 64 * 1024;
    std::vector<std::unique_ptr<Alignment>> alignments;
    while (parser->parse_objects(alignments, size_in_bytes)) {
    }

    EXPECT_EQ(48U, alignments.size());
}

TEST_F(BioparserSamTest, FormatError) {

    SetUp(bioparser_test_data_path + "sample.paf");

    std::vector<std::unique_ptr<Alignment>> alignments;

    EXPECT_DEATH(parser->parse_objects(alignments, -1),
        ".bioparser::SamParser. error: invalid file format!");
}
