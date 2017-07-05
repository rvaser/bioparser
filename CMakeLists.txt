cmake_minimum_required(VERSION 3.2)
project(bioparser)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

option(bioparser_build_tests "Build bioparser google tests" OFF)

add_library(bioparser INTERFACE)

target_include_directories(bioparser INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if (bioparser_build_tests)
    file(COPY test/data DESTINATION ${PROJECT_BINARY_DIR})
    set(bioparser_test_data_path ${PROJECT_BINARY_DIR}/data/)
    configure_file(
        "${PROJECT_SOURCE_DIR}/test/bioparser_test_config.h.in"
        "${PROJECT_BINARY_DIR}/config/bioparser_test_config.h"
    )
    include_directories(${PROJECT_BINARY_DIR}/config)

    add_executable(bioparser_test test/bioparser_test.cpp)

    add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)
    target_link_libraries(bioparser_test bioparser gtest_main)
endif(bioparser_build_tests)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/bioparser DESTINATION include)
