# bioparser

Bioparser is a c++ implementation of parsers for several bioinformatic formats. It consists of only one header file containing template parsers for FASTA, FASTQ, MHAP and PAF formats. Hopefully more formats will be added in the future.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+

## Usage

If you would like to add bioparser to your project, include the bioparser.hpp file while compiling and add -std=c++11 to your compiler flag list. For details on how to use the parsers, please look at the examples bellow:

    // define a class for sequences in FASTA format
    class ExampleClass {
        public:
            // required signature for the constructor
            ExampleClass(uint32_t object_id, const char* name, uint32_t name_length,
                const char* data, uint32_t data_length) {
                ...
            }
    };

    std::vector<std::unique_ptr<ExampleClass>> fasta_objects;
    auto fasta_reader = BIOPARSER::createReader<ExampleClass, BIOPARSER::FastaReader>(path_to_file);
    // read the whole file
    fasta_reader->read_objects(fasta_objects, -1);


    // define a class for sequences in FASTQ format
    class ExampleClass2 {
        public:
            // required signature for the constructor
            ExampleClass2(uint32_t object_id, const char* name, uint32_t name_length,
                const char* data, uint32_t data_length, const char* quality,
                uint32_t quality_length) {
                ...
            }
    };

    std::vector<std::unique_ptr<ExampleClass2>> fastq_objects;
    auto fastq_reader = BIOPARSER::createReader<ExampleClass2, BIOPARSER::FastqReader>(path_to_file2);
    // read a predefined size of bytes
    while (true) {
        auto status = fastq_reader->read_objects(fastq_objects, size_in_bytes);
        ...
        if (status == false) {
            break;
        }
    }


    // define a class for overlaps in MHAP format
    class ExampleClass3 {
        public:
            // required signature for the constructor
            ExampleClass3(uint32_t object_id, uint32_t a_id, uint32_t b_id, double eq_bases_perc, uint32_t minmers,
                uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
                uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length) {
                ...
            }
    };

    std::vector<std::unique_ptr<ExampleClass3>> mhap_objects;
    auto mhap_reader = BIOPARSER::createReader<ExampleClass3, BIOPARSER::MhapReader>(path_to_file3);
    mhap_reader->read_objects(mhap_objects, -1);


    // define a class for overlaps in PAF format or add a constructor to existing overlap class
    ...
        ExampleClass3(uint32_t object_id, const char* a_name, uint32_t a_name_length,
            uint32_t a_begin, uint32_t a_end, char orientation, const char* b_name,
            uint32_t b_name_length, uint32_t b_begin, uint32_t b_end,
            uint32_t matching_bases, uint32_t overlap_length,
            uint32_t quality {
            ...
        }
    ...

    std::vector<std::unique_ptr<ExampleClass3>> paf_objects;
    auto phap_reader = BIOPARSER::createReader<ExampleClass3, BIOPARSER::PafReader>(path_to_file4);
    paf_reader->read_objects(paf_objects, -1);

If your class has a **private** constructor with the required signature, format your classes in the following way:

    class ExampleClass {
        public:
            friend BIOPARSER::FastaReader<ExampleClass>;
        private:
            ExampleClass(...) {
                ...
            }
    }

    class ExampleClass2 {
        public:
            friend BIOPARSER::FastqReader<ExampleClass2>;
        private:
            ExampleClass2(...) {
                ...
            }
    }

    class ExampleClass3 {
        public:
            friend BIOPARSER::MhapReader<ExampleClass3>;
            friend BIOPARSER::PafReader<ExampleClass3>;
        private:
            ExampleClass3(...) {
                ...
            }
    }
