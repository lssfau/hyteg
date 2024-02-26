cmake_minimum_required (VERSION 3.20)

PROJECT ( hyteg )
enable_testing()
set_property(GLOBAL PROPERTY USE FOLDERS TRUE)

if(${CMAKE_VERSION} VERSION_GREATER "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${hyteg_SOURCE_DIR}/cmake/modules)
include(DetectPlatformFeatures)

option ( HYTEG_BUILD_WITH_MPFR        "Build with MPFR"                              OFF)
option ( HYTEG_BUILD_WITH_PETSC       "Build with PETSc"                             OFF)
option ( HYTEG_BUILD_WITH_TRILINOS    "Build with Trilinos"                          OFF)
option ( HYTEG_BUILD_WITH_ADIOS2      "Build with ADIOS2"                            OFF)
option ( HYTEG_BUILD_WITH_PYTHON3     "Build with PYTHON3"                           OFF)
option ( HYTEG_BUILD_WITH_AVX         "Use AVX intrinsics (in generated code)" ${HYTEG_PLATFORM_SUPPORTS_AVX})
option ( HYTEG_USE_GENERATED_KERNELS  "Use generated pystencils kernels if available" ON)
option ( HYTEG_GIT_SUBMODULE_AUTO     "Check submodules during build"                 ON)
option ( HYTEG_TERRANEO_MODULE        "Build TerraNeo module"                        OFF)
option ( HYTEG_DOWNLOAD_BOOST         "Download boost library; required for plates"  OFF)

set(WALBERLA_OPTIMIZE_FOR_LOCALHOST ON CACHE BOOL "Enable compiler optimizations spcific to localhost")
set(WALBERLA_LOG_SKIPPED            ON CACHE BOOL "Log skipped cmake targets"                         )
set(WALBERLA_BUILD_DOC              OFF )
set(WALBERLA_BUILD_WITH_BOOST       OFF )
set(WALBERLA_BUILD_TESTS            OFF )
set(WALBERLA_BUILD_BENCHMARKS       OFF )
set(WALBERLA_BUILD_TUTORIALS        OFF )
set(WALBERLA_BUILD_SHOWCASES        OFF )
set(WALBERLA_DIR ${hyteg_SOURCE_DIR}/walberla CACHE  PATH  "waLBerla path")

set(EIGEN_DIR ${hyteg_SOURCE_DIR}/eigen CACHE  PATH  "eigen path")
include_directories(SYSTEM ${EIGEN_DIR})
link_directories(${EIGEN_DIR})

include_directories(src)

unset(CMAKE_CXX_FLAGS_DEBUGOPTIMIZED CACHE)
find_package(waLBerla)

if ( HYTEG_BUILD_WITH_MPFR )
    # Extends cmake module path - so that FindMPFR.cmake in the current directory is found
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${EIGEN_DIR}/cmake)

    find_package(MPFR REQUIRED)
    add_library(MPFR::MPFR INTERFACE IMPORTED)
    set_target_properties(MPFR::MPFR PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDES}"
            INTERFACE_LINK_LIBRARIES "${MPFR_LIBRARIES}"
    )
    # add the path of the mpfr-wrapper to the include path
    include_directories(SYSTEM extern/mpfr-wrapper)
endif ()

if (NOT WALBERLA_DOUBLE_ACCURACY AND HYTEG_USE_GENERATED_KERNELS)
   set(HYTEG_USE_GENERATED_KERNELS OFF)
   message(WARNING "The generated constant stencil kernels are disabled when using WALBERLA_DOUBLE_ACCURACY=OFF")
endif ()

if (HYTEG_USE_GENERATED_KERNELS)
   message(STATUS "Using generated HyTeG kernels.")
else ()
   message(STATUS "Generated HyTeG kernels DISABLED! - Performance might not be optimal and some features might not be working correctly.")
endif ()

if(WALBERLA_CXX_COMPILER_IS_GNU)
    add_flag ( CMAKE_CXX_FLAGS "-Wno-comment" )
    add_flag ( CMAKE_CXX_FLAGS "-Wno-stringop-overflow" )
    if( NOT WALBERLA_DOUBLE_ACCURACY )
       add_flag ( CMAKE_CXX_FLAGS "-Wno-error=float-conversion" )
    endif()
endif()

include(FetchContent)
FetchContent_Declare(
      boost
      URL      https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.bz2
      URL_HASH SHA256=475d589d51a7f8b3ba2ba4eda022b170e562ca3b760ee922c146b6c65856ef39
)

if( HYTEG_DOWNLOAD_BOOST )
   FetchContent_Populate( boost )
   set(BOOST_ROOT ${hyteg_BINARY_DIR}/_deps/boost-src/ CACHE PATH "Path to boost root")
endif()

if(HYTEG_BUILD_WITH_PYTHON3)
   find_package (Python3 COMPONENTS Development)
endif()

# Add flags and libraries required to build terraneo submodule
if( HYTEG_TERRANEO_MODULE )

   find_package( Boost 1.48 )
   if( NOT Boost_FOUND )
      message( FATAL_ERROR "Boost is required for geometric calculations on the sphere but was not found.
      Please define BOOST_ROOT (e.g. cmake -DBOOST_ROOT=/local/boost/)
      or enable the cmake option HYTEG_DOWNLOAD_BOOST (e.g. cmake -DHYTEG_DOWNLOAD_BOOST=ON)")
   endif()

   # Unpack plate topology file
   if( NOT EXISTS ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson )
     if( EXISTS ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson.bz2 )
       message( STATUS "Unpacking plate topologies file" )
       execute_process(COMMAND bunzip2 "--keep" "topologies0-100Ma.geojson.bz2" WORKING_DIRECTORY ${hyteg_SOURCE_DIR}/data/terraneo/plates )
     else()
       message( WARNING "Missing plate topologies file: ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson.bzip2" )
     endif()
   endif()
endif()

# Clang complains about braces in correct init list
# see: https://stackoverflow.com/questions/31555584/why-is-clang-warning-suggest-braces-around-initialization-of-subobject-wmis
if(WALBERLA_CXX_COMPILER_IS_CLANG)
    add_flag( CMAKE_CXX_FLAGS "-Wno-missing-braces" )
    add_flag( CMAKE_CXX_FLAGS "-Wno-parentheses-equality" )
    if( NOT WALBERLA_DOUBLE_ACCURACY )
       add_flag ( CMAKE_CXX_FLAGS "-Wno-error=implicit-float-conversion" )
    endif()
endif()

# MSVC generates so many warnings that the CI log gets too large to complete the job successfully
# Therefore we suppress some warnings. This should be fixed at some point
if(WALBERLA_CXX_COMPILER_IS_MSVC)
    add_flag( CMAKE_CXX_FLAGS "/wd4100" )
    #suppress warning about missing case labels which are not present in generated kernels
    add_flag( CMAKE_CXX_FLAGS "/wd4065" )
endif()

if ( HYTEG_BUILD_WITH_PETSC )
    find_package(PETSc REQUIRED)
    add_library(PETSc::PETSc INTERFACE IMPORTED)
    set_target_properties(PETSc::PETSc PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDES}"
          INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES}"
          )
    set(WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE)
    message(STATUS "WALBERLA_BUILD_WITH_MPI was force set to ON")
    set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDES})
    include(CheckSymbolExists)
    unset(HYTEG_PETSC_BUILT_WITH_HDF5 CACHE)
    check_symbol_exists(PETSC_HAVE_HDF5 petscconf.h  HYTEG_PETSC_BUILT_WITH_HDF5)
endif()

if ( HYTEG_BUILD_WITH_TRILINOS )
   find_package( Trilinos REQUIRED )
   message(STATUS "Found Trilinos! Trilinos_DIR = ${Trilinos_DIR} ")
   INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
   list ( APPEND SERVICE_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARY_DIRS})
endif()

if ( HYTEG_BUILD_WITH_ADIOS2 )
  if ( WALBERLA_BUILD_WITH_MPI )
    find_package( ADIOS2 REQUIRED COMPONENTS CXX MPI )
  else()
    find_package( ADIOS2 REQUIRED COMPONENTS CXX )
  endif()
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

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/terraneo)
   file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/terraneo")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/terraneo/plates)
   file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/terraneo/plates")
endif()

waLBerla_link_files_to_builddir(data/meshes/*.msh)
waLBerla_link_files_to_builddir(data/meshes/3D/*.msh)
waLBerla_link_files_to_builddir(data/param/*.prm)
waLBerla_link_files_to_builddir(data/kerncraftMachineFiles/*.yml)
waLBerla_link_files_to_builddir(doc/images/*.png)
waLBerla_link_files_to_builddir(data/terraneo/*.json)
waLBerla_link_files_to_builddir(data/terraneo/plates/*.geojson)
waLBerla_link_files_to_builddir(data/terraneo/plates/*.rot)

add_subdirectory( src )
add_subdirectory( apps )
add_subdirectory( tests )
add_subdirectory( tutorials )


############################################################################################################################
# Documentation Generation
#
# Build documentation using Doxygen (www.doxygen.org)
############################################################################################################################
  find_package ( Doxygen )
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

  # The following checks all our "feature macros" and defines specific
  # variables for Doxygen, so that we can add the acitvated feature macros
  # to the list of PREDEFINED macros in Doxygen.
  #
  # Yepp, this looks akward. If you have any idea on how to solve this in
  # a more elegant fashion, please feel free to change this.
  if( WALBERLA_BUILD_WITH_MPI )
    set( DOXYGEN_WALBERLA_BUILD_WITH_MPI "WALBERLA_BUILD_WITH_MPI" )
  endif()

  if ( HYTEG_BUILD_WITH_ADIOS2 )
    set( DOXYGEN_HYTEG_BUILD_WITH_ADIOS2 "HYTEG_BUILD_WITH_ADIOS2" )
  endif()

  if ( HYTEG_BUILD_WITH_AVX )
    set( DOXYGEN_HYTEG_BUILD_WITH_AVX "HYTEG_BUILD_WITH_AVX" )
  endif()

  if ( HYTEG_BUILD_WITH_MPFR )
    set( DOXYGEN_HYTEG_BUILD_WITH_MPFR "HYTEG_BUILD_WITH_MPFR" )
  endif()

  if ( HYTEG_BUILD_WITH_PETSC )
    set( DOXYGEN_HYTEG_BUILD_WITH_PETSC "HYTEG_BUILD_WITH_PETSC" )
  endif()

  if ( HYTEG_BUILD_WITH_PYTHON3 )
    set( DOXYGEN_HYTEG_BUILD_WITH_PYTHON3 "HYTEG_BUILD_WITH_PYTHON3" )
  endif()

  if ( HYTEG_BUILD_WITH_TRILINOS )
    set( DOXYGEN_HYTEG_BUILD_WITH_TRILINOS "HYTEG_BUILD_WITH_TRILINOS" )
  endif()

  if ( HYTEG_USE_GENERATED_KERNELS )
    set( DOXYGEN_HYTEG_USE_GENERATED_KERNELS "HYTEG_USE_GENERATED_KERNELS" )
  endif()

  if ( HYTEG_TERRANEO_MODULE )
    set( DOXYGEN_HYTEG_TERRANEO_MODULE "HYTEG_TERRANEO_MODULE" )
  endif()

  configure_file ( ${hyteg_SOURCE_DIR}/doc/doxygen.config ${hyteg_BINARY_DIR}/doc/doxygen.cfg @ONLY )

  add_custom_target ( doc_hyteg   ${DOXYGEN_EXECUTABLE} ${hyteg_BINARY_DIR}/doc/doxygen.cfg
    COMMENT "Generating API documentation with Doxygen" VERBATIM )

else ()

  message( WARNING "Doxygen not found. You will not be able to build documentation." )

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
set(HYTEG_COMPILER_FLAGS "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} ${CMAKE_CXX_FLAGS}")
set(HYTEG_MPI_VERSION ${HYTEG_MPI_VERSION})
string(REPLACE "\"" "\\\"" HYTEG_COMPILER_FLAGS ${HYTEG_COMPILER_FLAGS})

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

message(STATUS "Detected Byte Order is: " ${CMAKE_CXX_BYTE_ORDER} )

############################################################################################################################
