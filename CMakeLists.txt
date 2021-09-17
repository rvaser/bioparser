cmake_minimum_required(VERSION 3.11)

project(bioparser VERSION 3.0.15
                  LANGUAGES CXX
                  DESCRIPTION "Bioparser is a c++ header only parsing library for several formats in bioinformatics (FASTA/Q, MHAP/PAF/SAM), with support for zlib compressed files.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(bioparser_main_project ON)
endif ()
option(bioparser_install "Generate install target" ${bioparser_main_project})
option(bioparser_build_tests "Build unit tests" ${bioparser_main_project})

find_package(ZLIB 1.2.8 REQUIRED)

if (bioparser_build_tests)
  find_package(biosoup 0.10.0 QUIET)
  find_package(GTest 1.10.0 QUIET)

  if (NOT biosoup_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      biosoup
      GIT_REPOSITORY https://github.com/rvaser/biosoup
      GIT_TAG 0.10.0)

    FetchContent_GetProperties(biosoup)
    if (NOT biosoup_POPULATED)
      FetchContent_Populate(biosoup)
      add_subdirectory(
        ${biosoup_SOURCE_DIR}
        ${biosoup_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()

  if (NOT GTest_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(bioparser INTERFACE)
add_library(${PROJECT_NAME}::bioparser ALIAS bioparser)

target_include_directories(bioparser INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(bioparser INTERFACE
  ZLIB::ZLIB)

if (bioparser_install)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS bioparser
    EXPORT ${PROJECT_NAME}Targets)
  install(
    DIRECTORY include/bioparser
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif ()

if (bioparser_build_tests)
  add_executable(bioparser_test
    test/parser_test.cpp
    test/fasta_parser_test.cpp
    test/fastq_parser_test.cpp
    test/mhap_parser_test.cpp
    test/paf_parser_test.cpp
    test/sam_parser_test.cpp)

  target_link_libraries(bioparser_test
    bioparser
    biosoup::biosoup
    GTest::Main)

  target_compile_definitions(bioparser_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")
endif ()
