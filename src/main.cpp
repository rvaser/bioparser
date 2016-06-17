#include <stdio.h>
#include <stdlib.h>

#include "bioparser.hpp"

using namespace BIOPARSER;

int main(int argc, char** argv) {

    auto reader = createReader<uint32_t>("Makefile", FileType::kFASTA);

    std::vector<std::unique_ptr<uint32_t>> test;
    reader->read_objects(test, 1000);

    return 0;
}
