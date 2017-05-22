cmake_minimum_required (VERSION 2.8)

PROJECT ( tinyhhg )
enable_testing()
option ( HHG_ENABLE_LIKWID                "Compile with LIKWID" )

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

include_directories ( extern/fmt-3.0.1 )
include_directories ( src )
#include_directories ( src/tinyhhg_core )
#include_directories ( ${MPI_CXX_INCLUDE_PATH})
#include_directories ( ${tinyhhg_BINARY_DIR}/src/tinyhhg_core )

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${tinyhhg_SOURCE_DIR} )
find_package( waLBerla )


add_subdirectory (extern/fmt-3.0.1)
add_subdirectory( src )
add_subdirectory( apps )
add_subdirectory( tests )