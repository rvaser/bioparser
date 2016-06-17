#include <stdio.h>
#include <stdlib.h>

#include "bioparser.hpp"

using namespace BIOPARSER;

class Read {
public:
    Read(const std::string& name, const std::string& data)
            : name_(name), data_(data){
    }
    std::string name_;
    std::string data_;
};

int main(int argc, char** argv) {

    auto reader = createReader<Read>(argv[1], FileType::kFASTA);

    std::vector<std::unique_ptr<Read>> test;
    reader->read_objects(test, 10000000000000);

    return 0;
}
