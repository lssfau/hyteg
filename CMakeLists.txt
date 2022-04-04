cmake_minimum_required (VERSION 3.14)

PROJECT ( hyteg )
enable_testing()
set_property(GLOBAL PROPERTY USE FOLDERS TRUE)

option ( HYTEG_BUILD_WITH_PETSC       "Build with PETSc"                             OFF)
option ( HYTEG_BUILD_WITH_TRILINOS    "Build with Trilinos"                          OFF)
option ( HYTEG_USE_GENERATED_KERNELS  "Use generated pystencils kernels if available" ON)
option ( HYTEG_GIT_SUBMODULE_AUTO     "Check submodules during build"                 ON)

set(WALBERLA_OPTIMIZE_FOR_LOCALHOST ON  CACHE BOOL "Enable compiler optimizations spcific to localhost")
set(WALBERLA_LOG_SKIPPED            ON CACHE BOOL "Log skipped cmake targets"                         )
set(WALBERLA_DIR ${hyteg_SOURCE_DIR}/walberla CACHE  PATH  "waLBerla path")

set(EIGEN_DIR ${hyteg_SOURCE_DIR}/eigen CACHE  PATH  "eigen path")
include_directories ( ${EIGEN_DIR} )
link_directories    ( ${EIGEN_DIR} )

include_directories ( src )

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${hyteg_SOURCE_DIR}/cmake/modules )

if ( HYTEG_USE_GENERATED_KERNELS )
    message(STATUS "Using generated HyTeG kernels.")
else()
    message(STATUS "Generated HyTeG kernels DISABLED! - Performance might not be optimal and some features might not be working correctly.")
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

# MSVC generates so many warnings that the CI log gets too large to complete the job successfully
# Therefore we suppress some warnings. This should be fixed at some point
if(WALBERLA_CXX_COMPILER_IS_MSVC)
    add_flag( CMAKE_CXX_FLAGS "/wd4100" )
    #suppress warning about missing case labels which are not present in generated kernels
    add_flag( CMAKE_CXX_FLAGS "/wd4065" )
endif()

if ( HYTEG_BUILD_WITH_PETSC )
    find_package( PETSc REQUIRED )
    include_directories ( ${PETSC_INCLUDES} )
    link_directories    ( ${PETSC_LIBRARIES} )
    list ( APPEND SERVICE_LIBS ${PETSC_LIBRARIES} )
    set(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE)
    message(STATUS "WALBERLA_BUILD_WITH_MPI was force set to ON")
endif()

if ( HYTEG_BUILD_WITH_TRILINOS)
   find_package( Trilinos REQUIRED )
   message(STATUS "Found Trilinos! Trilinos_DIR = ${Trilinos_DIR} ")
   INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
   list ( APPEND SERVICE_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARY_DIRS})
endif()


configure_file ( ${hyteg_SOURCE_DIR}/src/hyteg/HytegDefinitions.in.hpp
                 src/hyteg/HytegDefinitions.hpp )

include_directories( ${hyteg_BINARY_DIR}/src )


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

if(NOT EXISTS ${CMAKE_BINARY_DIR}/doc/images)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/doc/images")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/kerncraftMachineFiles)
   file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/kerncraftMachineFiles")
endif()

waLBerla_link_files_to_builddir(data/meshes/*.msh)
waLBerla_link_files_to_builddir(data/meshes/3D/*.msh)
waLBerla_link_files_to_builddir(data/param/*.prm)
waLBerla_link_files_to_builddir(data/kerncraftMachineFiles/*.yml)
waLBerla_link_files_to_builddir(doc/images/*.png)

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
   set ( DOXYGEN_HTML_HEADER ${hyteg_SOURCE_DIR}/doc/header.html )
   set ( DOXYGEN_HTML_FOOTER ${hyteg_SOURCE_DIR}/doc/footer.html )
   set ( DOXYGEN_HTML_OUTPUT "html" )

   configure_file ( ${hyteg_SOURCE_DIR}/doc/doxygen.config ${hyteg_BINARY_DIR}/doc/doxygen.cfg @ONLY )

   add_custom_target ( doc_hyteg   ${DOXYGEN_EXECUTABLE} ${hyteg_BINARY_DIR}/doc/doxygen.cfg
                       COMMENT "Generating API documentation with Doxygen" VERBATIM )

endif ( )
############################################################################################################################


############################################################################################################################
# Git commit hash and submodules
############################################################################################################################

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
   # Update submodules as needed

   if(GIT_SUBMODULE)
      message(STATUS "Submodule update")
      execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()
   endif()
   # Get the current working branch
   execute_process(
         COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
         OUTPUT_VARIABLE GIT_BRANCH
         OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   # Get the latest abbreviated commit hash of the working branch
   execute_process(
         COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
         OUTPUT_VARIABLE GIT_COMMIT_HASH
         OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
   add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")

   configure_file(
         ${CMAKE_SOURCE_DIR}/src/hyteg/Git.in.hpp
         ${CMAKE_BINARY_DIR}/src/hyteg/Git.hpp
   )
else()
   message(FATAL_ERROR "git executalbe could not be found")
endif()



############################################################################################################################

############################################################################################################################
# Build info
############################################################################################################################

string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )

set      ( HYTEG_BUILD_TYPE      ${CMAKE_BUILD_TYPE} )
set      ( HYTEG_COMPILER_INFO   "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}" )
set      ( HYTEG_COMPILER_FLAGS  "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${CMAKE_CXX_FLAGS}" )
string   ( REPLACE "\"" "\\\"" HYTEG_COMPILER_FLAGS ${HYTEG_COMPILER_FLAGS} )

configure_file(
        ${CMAKE_SOURCE_DIR}/src/hyteg/BuildInfo.in.hpp
        ${CMAKE_BINARY_DIR}/src/hyteg/BuildInfo.hpp
)


# this returns all targets within one folder
function(get_all_targets _result _dir)
   get_property(_subdirs DIRECTORY "${_dir}" PROPERTY SUBDIRECTORIES)
   foreach (_subdir IN LISTS _subdirs)
      get_all_targets(${_result} "${_subdir}")
   endforeach ()
   get_property(_sub_targets DIRECTORY "${_dir}" PROPERTY BUILDSYSTEM_TARGETS)
   set(${_result} ${${_result}} ${_sub_targets} PARENT_SCOPE)
endfunction()

get_all_targets(_tests_targets "tests")
add_custom_target(hyteg_all_tests)
add_dependencies(hyteg_all_tests ${_tests_targets})

get_all_targets(_apps_targets "apps")
get_all_targets(_tutorials_targets "tutorials")
add_custom_target(hyteg_all_apps)
add_dependencies(hyteg_all_apps ${_apps_targets} ${_tutorials_targets})

add_custom_target(hyteg_all)
add_dependencies(hyteg_all ${_apps_targets} ${_tutorials_targets} ${_tests_targets})
############################################################################################################################
