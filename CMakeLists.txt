cmake_minimum_required(VERSION 3.16.3)
project(pufparser)

set(CMAKE_CXX_STANDARD 17)

add_executable(pufparser main.cpp DataParser.cpp DataParser.h)
