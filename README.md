# Bioparser

![image](https://travis-ci.org/rvaser/bioparser.svg?branch=master)

Bioparser is a c++ implementation of parsers for several bioinformatics formats. It consists of only one header file containing template parsers for FASTA, FASTQ, MHAP and PAF formats.

## Dependencies

### Linux

Application uses following software:

1. gcc 4.8+ or clang 3.4+
2. (optional) cmake 3.2+

## Usage

If you would like to add bioparser to your project, add `-Iinclude/` and `-std=c++11` while compiling and include `bioparser/bioparser.hpp` in your desired source files. Alternatively, add the project to your CMakeLists.txt file with the `add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)` and `target_link_libraries(your_exe bioparser)` commands.

For details on how to use the parsers in your code, please look at the examples bellow:

```cpp
// define a class for sequences in FASTA format
class ExampleClass {
    public:
        // required signature for the constructor
        ExampleClass(uint64_t object_id,
            const char* name,
            uint32_t name_length,
            const char* data,
            uint32_t data_length) {
            // your implementation
        }
};

std::vector<std::unique_ptr<ExampleClass>> fasta_objects;
auto fasta_reader = bioparser::createReader<ExampleClass, bioparser::FastaReader>(path_to_file);
// read the whole file
fasta_reader->read_objects(fasta_objects, -1);

// define a class for sequences in FASTQ format
class ExampleClass2 {
    public:
        // required signature for the constructor
        ExampleClass2(uint64_t object_id,
            const char* name,
            uint32_t name_length,
            const char* data,
            uint32_t data_length,
            const char* quality,
            uint32_t quality_length) {
            // your implementation
        }
};

std::vector<std::unique_ptr<ExampleClass2>> fastq_objects;
auto fastq_reader = bioparser::createReader<ExampleClass2, bioparser::FastqReader>(path_to_file2);
// read a predefined size of bytes
uint64_t size_in_bytes = 500 * 1024 * 1024; // 500 MB
while (true) {
    auto status = fastq_reader->read_objects(fastq_objects, size_in_bytes);
    // do some work with objects
    if (status == false) {
        break;
    }
}

// define a class for overlaps in MHAP format
class ExampleClass3 {
    public:
        // required signature for the constructor
        ExampleClass3(uint64_t object_id,
            uint32_t a_id,
            uint32_t b_id,
            double eq_bases_perc,
            uint32_t minmers,
            uint32_t a_rc,
            uint32_t a_begin,
            uint32_t a_end,
            uint32_t a_length,
            uint32_t b_rc,
            uint32_t b_begin,
            uint32_t b_end,
            uint32_t b_length) {
            // your implementation
        }
};

std::vector<std::unique_ptr<ExampleClass3>> mhap_objects;
auto mhap_reader = bioparser::createReader<ExampleClass3, bioparser::MhapReader>(path_to_file3);
mhap_reader->read_objects(mhap_objects, -1);

// define a class for overlaps in PAF format or add a constructor to existing overlap class
ExampleClass3::ExampleClass3(uint64_t object_id,
    const char* a_name,
    uint32_t a_name_length,
    uint32_t a_begin,
    uint32_t a_end,
    char orientation,
    const char* b_name,
    uint32_t b_name_length,
    uint32_t b_begin,
    uint32_t b_end,
    uint32_t matching_bases,
    uint32_t overlap_length,
    uint32_t quality {
    // your implementation
}

std::vector<std::unique_ptr<ExampleClass3>> paf_objects;
auto phap_reader = bioparser::createReader<ExampleClass3, bioparser::PafReader>(path_to_file4);
paf_reader->read_objects(paf_objects, -1);
```
If your class has a **private** constructor with the required signature, format your classes in the following way:

```cpp
class ExampleClass {
    public:
        friend bioparser::FastaReader<ExampleClass>;
    private:
        ExampleClass(...) {
            ...
        }
}

class ExampleClass2 {
    public:
        friend bioparser::FastqReader<ExampleClass2>;
    private:
        ExampleClass2(...) {
            ...
        }
}

class ExampleClass3 {
    public:
        friend bioparser::MhapReader<ExampleClass3>;
        friend bioparser::PafReader<ExampleClass3>;
    private:
        ExampleClass3(...) {
            ...
        }
}
```
