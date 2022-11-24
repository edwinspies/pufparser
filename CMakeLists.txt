cmake_minimum_required(VERSION 3.24)
project(pufparser)

set(CMAKE_CXX_STANDARD 17)

add_executable(pufparser main.cpp DataParser.cpp DataParser.h Google_tests/lib/test_DataParser.cpp)
