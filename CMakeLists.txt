cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

set (CMAKE_CXX_STANDARD 11)
# set (CMAKE_CXX_FLAGS "-O2")
set (CMAKE_CXX_FLAGS "-g -O0")

add_subdirectory (extern/fmt-3.0.1)
include_directories ("${PROJECT_SOURCE_DIR}/extern/fmt-3.0.1")

include_directories ("${PROJECT_SOURCE_DIR}/include")
add_subdirectory (src)

add_executable (tinytest test/test.cpp)
target_link_libraries (tinytest tinyhhg fmt ${MPI_LIBRARIES})

add_executable (tinytest_cg test/cg_test.cpp)
target_link_libraries (tinytest_cg tinyhhg fmt ${MPI_LIBRARIES})

add_executable (tinytest_fmg test/fmg_test.cpp)
target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_LIBRARIES})