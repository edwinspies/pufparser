cmake_minimum_required(VERSION 3.7.2)
project(pufparser)

set(CMAKE_CXX_STANDARD 11)

add_executable(pufparser main.cpp DataParser.cpp DataParser.h)
