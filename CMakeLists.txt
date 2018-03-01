cmake_minimum_required (VERSION 2.8)

PROJECT ( tinyhhg )
enable_testing()

option ( HHG_BUILD_WITH_PETSC             "Build with PETSc"                       OFF)
option ( HHG_REMOVE_WALBERLA              "Remove unnecessary walberla components" OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include_directories ( src )

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${tinyhhg_SOURCE_DIR}/cmake/modules )

if ( HHG_BUILD_WITH_PETSC )
    find_package( PETSc )
    include_directories ( ${PETSC_INCLUDE_DIRS} )
    link_directories    ( ${PETSC_LIBRARY_DIR} )
    list ( APPEND SERVICE_LIBS ${PETSC_LIBRARIES} )
    set(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE)
    message(STATUS "WALBERLA_BUILD_WITH_MPI was force set to ON")
endif()

if( HHG_REMOVE_WALBERLA )
    execute_process(COMMAND ${tinyhhg_SOURCE_DIR}/cmake/patchWalberla.sh WORKING_DIRECTORY ${WALBERLA_DIR})
endif()

find_package( waLBerla )

# We need to add the c++14 flag since some older cmake versions do not work with intel and CXX_STANDARD
if(WALBERLA_CXX_COMPILER_IS_INTEL)
    add_flag ( CMAKE_CXX_FLAGS "-std=c++14" )
endif()

# Clang complains about braces in correct init list
# see: https://stackoverflow.com/questions/31555584/why-is-clang-warning-suggest-braces-around-initialization-of-subobject-wmis
if(WALBERLA_CXX_COMPILER_IS_CLANG)
    add_flag( CMAKE_CXX_FLAGS "-Wno-missing-braces" )
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

waLBerla_link_files_to_builddir(data/meshes/*.msh)
waLBerla_link_files_to_builddir(data/meshes/3D/*.msh)
waLBerla_link_files_to_builddir(data/param/*.prm)

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
