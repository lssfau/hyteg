cmake_minimum_required (VERSION 2.8)

PROJECT ( tinyhhg )
enable_testing()

option ( HHG_BUILD_WITH_PETSC       "Build with PETSc"                             OFF)
option ( HHG_BUILD_WITH_EIGEN       "Build with Eigen"                             OFF)
option ( HHG_USE_GENERATED_KERNELS  "Use generated pystencils kernels if available" ON)
option ( HHG_P1_COLORING            "Use multicolored P1 layout EXPERIMENTAL (tets only) - only toggle when using only P1 elements" OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(WALBERLA_OPTIMIZE_FOR_LOCALHOST ON CACHE BOOL "Enable compiler optimizations spcific to localhost")
set(WALBERLA_DIR ${tinyhhg_SOURCE_DIR}/walberla CACHE  PATH  "waLBerla path")


include_directories ( src )

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${tinyhhg_SOURCE_DIR}/cmake/modules )

if ( HHG_P1_COLORING )
    message(STATUS "Using multicolor layout for P1 elements! Do not use other discretizations at all.")
endif()

if ( HHG_USE_GENERATED_KERNELS )
    message(STATUS "Using generated HyTeG kernels.")
else()
    message(STATUS "Generated HyTeG kernels DISABLED! - Performance might not be optimal and some features might not be working correctly.")
endif()

if ( HHG_BUILD_WITH_PETSC )
    find_package( PETSc )
    include_directories ( ${PETSC_INCLUDES} )
    link_directories    ( ${PETSC_LIBRARIES} )
    list ( APPEND SERVICE_LIBS ${PETSC_LIBRARIES} )
    set(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE)
    message(STATUS "WALBERLA_BUILD_WITH_MPI was force set to ON")
endif()

if ( HHG_BUILD_WITH_EIGEN )
    find_package( Eigen3 REQUIRED )
    include_directories ( ${EIGEN3_INCLUDE_DIRS} )
    link_directories    ( ${EIGEN3_LIBRARY_DIR} )
endif()

find_package( waLBerla )

if(WALBERLA_CXX_COMPILER_IS_GNU)
    add_flag ( CMAKE_CXX_FLAGS "-Wno-comment" )
    add_flag ( CMAKE_CXX_FLAGS "-Wno-stringop-overflow" )
endif()

# We need to add the c++14 flag since some older cmake versions do not work with intel and CXX_STANDARD
if(WALBERLA_CXX_COMPILER_IS_INTEL)
    add_flag ( CMAKE_CXX_FLAGS "-std=c++14" )
endif()

# Clang complains about braces in correct init list
# see: https://stackoverflow.com/questions/31555584/why-is-clang-warning-suggest-braces-around-initialization-of-subobject-wmis
if(WALBERLA_CXX_COMPILER_IS_CLANG)
    add_flag( CMAKE_CXX_FLAGS "-Wno-missing-braces" )
    add_flag( CMAKE_CXX_FLAGS "-Wno-parentheses-equality" )
endif()

# MSVC generates so many warnings that the CI log gets too large to complete the job successfully - therefore we suppress a prominent one.
if(WALBERLA_CXX_COMPILER_IS_MSVC)
    add_flag( CMAKE_CXX_FLAGS "/wd4100" )
endif()

configure_file ( ${tinyhhg_SOURCE_DIR}/src/tinyhhg_core/HHGDefinitions.in.hpp
                 src/tinyhhg_core/HHGDefinitions.hpp )

include_directories( ${tinyhhg_BINARY_DIR}/src )


if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/meshes)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/meshes")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/meshes/3D)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/meshes/3D")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/param)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/param")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/output")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/doc/ExtraPics)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/doc/ExtraPics")
endif()

waLBerla_link_files_to_builddir(data/meshes/*.msh)
waLBerla_link_files_to_builddir(data/meshes/3D/*.msh)
waLBerla_link_files_to_builddir(data/param/*.prm)
waLBerla_link_files_to_builddir(doc/ExtraPics/*.png)

add_subdirectory( src )
add_subdirectory( apps )
add_subdirectory( tests )
add_subdirectory( tutorials )


############################################################################################################################
# Documentation Generation
#
# Build documentation using Doxygen (www.doxygen.org)
############################################################################################################################
find_package ( Doxygen  )
find_package ( HTMLHelp )

if ( HTML_HELP_COMPILER EQUAL "" )
   set ( HTML_HELP_FOUND "NO" )
else ( )
   set ( HTML_HELP_FOUND "YES" )
endif ( )

if ( DOXYGEN_FOUND )
   set ( DOXYGEN_HTML_HEADER ${tinyhhg_SOURCE_DIR}/doc/header.html )
   set ( DOXYGEN_HTML_FOOTER ${tinyhhg_SOURCE_DIR}/doc/footer.html )
   set ( DOXYGEN_HTML_OUTPUT "html" )

   configure_file ( ${tinyhhg_SOURCE_DIR}/doc/doxygen.config ${tinyhhg_BINARY_DIR}/doc/doxygen.cfg @ONLY )

   add_custom_target ( doc_tinyhhg   ${DOXYGEN_EXECUTABLE} ${tinyhhg_BINARY_DIR}/doc/doxygen.cfg
                       COMMENT "Generating API documentation with Doxygen" VERBATIM )

endif ( )
############################################################################################################################
