#include <stdio.h>
#include <stdlib.h>

#include "bioparser.hpp"

using namespace BIOPARSER;

class Read {
public:
    Read(const char* name, uint32_t name_length, const char* data, uint32_t data_length)
            : name_(name, name_length), data_(data, data_length), quality_() {
    }

    Read(const char* name, uint32_t name_length, const char* data, uint32_t data_length,
        const char* quality, uint32_t quality_length)
            : name_(name, name_length), data_(data, data_length),
            quality_(quality, quality_length) {
    }

    std::string name_;
    std::string data_;
    std::string quality_;
};

int main(int argc, char** argv) {

    auto reader = createReader<Read>(argv[1], FileType::kFASTQ);

    std::vector<std::unique_ptr<Read>> test;
    reader->read_objects(test, 10000000000000);

    return 0;
}
