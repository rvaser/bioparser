# Bioparser

[![Latest GitHub release](https://img.shields.io/github/release/rvaser/bioparser.svg)](https://github.com/rvaser/bioparser/releases/latest)
![Build status for gcc/clang](https://github.com/rvaser/bioparser/actions/workflows/bioparser.yml/badge.svg)

Bioparser is a c++ header only parsing library for several bioinformatics formats (FASTA/Q, MHAP/PAF/SAM), with support for zlib compressed files.

## Usage

To build bioparser run the following commands:
```bash
git clone https://github.com/rvaser/bioparser && cd bioparser && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```
which will create install targets and unit tests. Running `make install` will create a package on your system that can be searched and linked with:
```cmake
find_package(bioparser)
target_link_libraries(<target> bioparser::bioparser)
```
On the other hand, you can include bioparser as a submodule and add it to your project with the following:
```cmake
if (NOT TARGET bioparser)
  add_subdirectory(<path_to_submodules>/bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<target> bioparser::bioparser)
```

If you are not using CMake, include the appropriate header file directly to your project and link with zlib.

#### Build options

- `bioparser_install`: generate install target
- `bioparser_build_tests`: build unit tests

#### Dependencies
- gcc 4.8+ | clang 3.5+
- zlib 1.2.8+
- (optional) cmake 3.11+

###### Hidden
- (bioparser_test) rvaser/biosoup 0.10.0
- (bioparser_test) google/googletest 1.10.0

## Examples

#### FASTA parser

```cpp
#include "bioparser/fasta_parser.hpp"

struct Sequence {  // or any other name
 public:
  Sequence(  // required arguments
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(path);

// parse whole file
auto s = p->Parse(-1);
```

#### FASTQ parser

```cpp
#include "bioparser/fastq_parser.hpp"

struct Sequence {  // or any other name
 public:
  Sequence(  // required arguments
      const char*, std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(path);

// parse in chunks
std::vector<std::unique_ptr<Sequence>> s;
while (true) {
  auto c = p->Parse(1ULL << 30);  // 1 GB
  if (c.empty()) {
    break;
  }
  s.insert(
      s.end(),
      std::make_move_iterator(c.begin()),
      std::make_move_iterator(c.end()));
}
```

#### MHAP parser

```cpp
#include "bioparser/mhap_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      std::uint64_t,
      std::uint64_t,
      double error,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::MhapParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

#### PAF parser

```cpp
#include "bioparser/paf_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      char,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::PafParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

#### SAM parser

```cpp
#include "bioparser/sam_parser.hpp"

struct Overlap {  // or any other name
 public:
  Overlap(  // required arguments
      const char*, std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      const char*, std::uint32_t,
      const char*, std::uint32_t) {
    // implementation
  }
}
auto p = bioparser::Parser<Overlap>::Create<bioparser::SamParser>(path);

// parse whole file
auto o = p->Parse(-1);
```

**Note**: If your class has a private constructor, add one of the following lines to your class definition:

```cpp
friend bioparser::FastaParser<Sequence>;
friend bioparser::FastqParser<Sequence>;
friend bioparser::MhapParser<Overlap>;
friend bioparser::PafParser<Overlap>;
friend bioparser::SamParser<Overlap>;
```

## Acknowledgement

This work has been supported in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
