# bioparser

Bioparser is a c++ implementation of a parser for bioinformatic formats. It consists of only one header file containing template parsers for FASTA, FASTQ and MHAP formats. Hopefully more formats will be added in the future.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+

## Usage

If you would like to add bioparser to your project, include the bioparser.hpp file while compiling and add -std=c++11 to your compiler flag list. All parser return a vector containing std::unique_ptr<YourClass>. For more details please look at the examples bellow:

    // user defined class for sequences in FASTA format
    class ExampleClass {
        public:
            ExampleClass(uint32_t object_id, const char* name, uint32_t name_length,
                const char* data, uint32_t data_length) {
                ...
            }
    };

    std::vector<std::unique_ptr<ExampleClass>> fasta_objects;
    auto fasta_reader = BIOPARSER::createFastaReader<ExampleClass>(path_to_file);
    // read a predefined size of bytes (size_in_bytes) or whole file (-1)
    fasta_reader->read_objects(fasta_objects, size_in_bytes);

    // user defined class for sequences in FASTQ format
    class ExampleClass2 {
        public:
            ExampleClass2(uint32_t object_id, const char* name, uint32_t name_length,
                const char* data, uint32_t data_length, const char* quality,
                uint32_t quality_length) {
                ...
            }
    };

    std::vector<std::unique_ptr<ExampleClass2>> fastq_objects;
    auto fastq_reader = BIOPARSER::createFastaReader<ExampleClass2>(path_to_file2);
    // read a predefined size of bytes (size_in_bytes) or whole file (-1)
    fastq_reader->read_objects(fastq_objects, size_in_bytes);


    // user defined class for overlaps in MHAP format
    class ExampleClass3 {
        ExampleClass3(uint32_t object_id, uint32_t a_id, uint32_t b_id,
            double eq_bases_perc, uint32_t minmers,
            uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
            uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length) {
            ...
        }
    };

    std::vector<std::unique_ptr<ExampleClass3>> mhap_objects;
    auto mhap_reader = BIOPARSER::createMhapReader<ExampleClass3>(path_to_file3);
    // read a predefined size of bytes (size_in_bytes) or whole file (-1)
    mhap_reader->read_objects(mhap_objects, size_in_bytes);

If your class has a private constructor with the required signature, add the following lines to your classes:

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
            friend BIOPARSEr::MhapReader<ExampleClass3>;
        private:
            ExampleClass3(...) {
                ...
            }
    }
