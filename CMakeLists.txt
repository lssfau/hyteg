cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})

set (CMAKE_CXX_STANDARD 11)

option ( HHG_ENABLE_LIKWID                "Compile with LIKWID" )

if ( NOT CMAKE_BUILD_TYPE )
  set ( CMAKE_BUILD_TYPE Release )
endif()

if ( HHG_ENABLE_LIKWID )
    find_library( LIKWID_LIB likwid HINTS $ENV{LIKWID_ROOT}/lib )
    find_path( LIKWID_INCLUDE_DIR likwid.h HINTS $ENV{LIKWID_ROOT}/include )
    message("HHG_ENABLE_LIKWID set")
    if ( LIKWID_LIB AND LIKWID_INCLUDE_DIR)
        include_directories( ${LIKWID_INCLUDE_DIR})
        add_definitions ( "-DLIKWID_PERFMON" )
    else()
        message(WARNING "likwid marker library not found. Set environment variable LIKWID_ROOT")
    endif()
endif()

add_subdirectory (extern/fmt-3.0.1)
include_directories ("${PROJECT_SOURCE_DIR}/extern/fmt-3.0.1")

include_directories ("${PROJECT_SOURCE_DIR}/include" ${MPI_CXX_INCLUDE_PATH})
add_subdirectory (src)

add_executable (tinytest test/test.cpp)
target_link_libraries (tinytest tinyhhg fmt ${MPI_CXX_LIBRARIES})

add_executable (tinytest_cg test/cg_test.cpp)
target_link_libraries (tinytest_cg tinyhhg fmt ${MPI_CXX_LIBRARIES})

add_executable (tinytest_fmg test/fmg_test.cpp)
target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})

add_executable (tinytest_chorin test/chorin_test.cpp)
target_link_libraries (tinytest_chorin tinyhhg fmt ${MPI_CXX_LIBRARIES})
