cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
project (tinyhhg)

set (CMAKE_CXX_STANDARD 11)

include_directories ("${PROJECT_SOURCE_DIR}/../include")

file(GLOB SOURCES *.cpp primitives/*.cpp)

add_library(tinyhhg ${SOURCES})