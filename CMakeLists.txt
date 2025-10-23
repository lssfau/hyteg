cmake_minimum_required (VERSION 3.20)

PROJECT ( hyteg )
enable_testing()
set_property(GLOBAL PROPERTY USE FOLDERS TRUE)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(${CMAKE_VERSION} VERSION_GREATER "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

if(${CMAKE_VERSION} VERSION_GREATER "3.27.0")
   cmake_policy(SET CMP0144 NEW)
endif()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.30.0")
   cmake_policy(SET CMP0167 OLD)
endif()

# Extends cmake module path - so that FindwaLBerla.cmake in the current directory is found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${hyteg_SOURCE_DIR}/cmake/modules)
include(DetectPlatformFeatures)

option ( HYTEG_BUILD_WITH_MPI            "Build with MPI"                                                             ON )
option ( HYTEG_BUILD_WITH_OPENMP         "Build with OpenMP"                                                         OFF )
option ( HYTEG_BUILD_WITH_PETSC          "Build with PETSc"                                                          OFF )
option ( HYTEG_BUILD_WITH_MPFR           "Build with MPFR"                                                           OFF )
option ( HYTEG_BUILD_WITH_TRILINOS       "Build with Trilinos"                                                       OFF )
option ( HYTEG_BUILD_WITH_ADIOS2         "Build with ADIOS2"                                                         OFF )
option ( HYTEG_BUILD_WITH_PYTHON3        "Build with PYTHON3"                                                        OFF )
option ( HYTEG_BUILD_WITH_AVX            "Use AVX/AVX2 intrinsics (in generated code)"    ${HYTEG_PLATFORM_SUPPORTS_AVX} )
option ( HYTEG_USE_GENERATED_KERNELS     "Use generated pystencils kernels if available"                              ON )
option ( HYTEG_GIT_SUBMODULE_AUTO        "Check submodules during build"                                              ON )
option ( HYTEG_TERRANEO_MODULE           "Build TerraNeo module"                                                     OFF )
option ( HYTEG_MANTLECONVECTION_APP      "Build Mantle convection app"                                               OFF )
option ( HYTEG_DONT_BUILD_DURING_CI      "Can be used to exclude build targets during CI"                             ON )
option ( HYTEG_DOWNLOAD_BOOST            "Download boost library; required for plates"                               OFF )
option ( HYTEG_USE_SIGNED_INT_FOR_ADIOS2 "Use int64_t when writing BP files for ParaView"                             ON )
# option ( HYTEG_OPTIMIZE_FOR_LOCALHOST    "Enable compiler optimizations specific to localhost"                        ON )

# Allow setting this value via CMake. Default should be 0, which is Doxygen's default, too.
set ( DOXYGEN_LOOKUP_CACHE_SIZE 0 CACHE STRING "See Doxygen documentation LOOKUP_CACHE_SIZE" )
mark_as_advanced( DOXYGEN_LOOKUP_CACHE_SIZE )

if ( HYTEG_BUILD_WITH_MPI )
   set(WALBERLA_BUILD_WITH_MPI ON)
else ()
   set( WALBERLA_BUILD_WITH_MPI OFF )
endif ()

set(WALBERLA_LOG_SKIPPED            ON CACHE BOOL "Log skipped cmake targets"                         )
set(WALBERLA_BUILD_DOC              OFF )
set(WALBERLA_BUILD_WITH_BOOST       OFF )
set(WALBERLA_BUILD_WITH_OPENMESH    OFF )
set(WALBERLA_BUILD_TESTS            OFF )
set(WALBERLA_BUILD_BENCHMARKS       OFF )
set(WALBERLA_BUILD_TUTORIALS        OFF )
set(WALBERLA_BUILD_SHOWCASES        OFF )
set(WALBERLA_DEPS_ERROR             ON  )
set(WALBERLA_DIR ${hyteg_SOURCE_DIR}/walberla CACHE  PATH  "waLBerla path")

set(EIGEN_DIR ${hyteg_SOURCE_DIR}/eigen CACHE  PATH  "eigen path")
include_directories(SYSTEM ${EIGEN_DIR})
link_directories(${EIGEN_DIR})

add_subdirectory( walberla EXCLUDE_FROM_ALL SYSTEM )
list( APPEND WALBERLA_MODULE_DIRS "${hyteg_SOURCE_DIR}/src" "${hyteg_SOURCE_DIR}/tests" )

# mark walberla include as system to prevent warning from walberla code
target_include_directories(core SYSTEM PUBLIC ${WALBERLA_DIR}/src)
target_include_directories( core SYSTEM PUBLIC ${hyteg_BINARY_DIR}/walberla/src )

add_flag ( CMAKE_CXX_FLAGS "-Wall -Wconversion -Wshadow" )
if ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" )
   add_flag( CMAKE_CXX_FLAGS "-Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-float-conversion -Wno-unused-function -Wno-unused-variable" )
   add_flag( CMAKE_CXX_FLAGS "-Wno-implicit-float-conversion -Wno-unused-local-typedef -Wno-delete-non-abstract-non-virtual-dtor" )
   add_flag( CMAKE_CXX_FLAGS "-Wno-unused-but-set-variable -Wno-implicit-int-conversion" )
elseif ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
   add_flag( CMAKE_CXX_FLAGS "-Wno-conversion -Wno-float-conversion -Wno-sign-compare -Wno-sign-conversion -Wno-unused-function" )
   add_flag( CMAKE_CXX_FLAGS "-Wno-unused-variable -Wno-unused-local-typedefs -Wno-unused-but-set-variable -Wno-comment" )
endif ()

# Treat warnings as errors
if ( WARNING_ERROR )
   if( WALBERLA_CXX_COMPILER_IS_GNU OR WALBERLA_CXX_COMPILER_IS_INTEL OR WALBERLA_CXX_COMPILER_IS_CLANG )
      add_flag ( CMAKE_CXX_FLAGS "-pedantic-errors -Werror" )
   elseif( WALBERLA_CXX_COMPILER_IS_MSVC )
      add_flag ( CMAKE_CXX_FLAGS "/WX" )
   elseif ( WALBERLA_CXX_COMPILER_IS_CRAY )
      add_flag ( CMAKE_CXX_FLAGS "-h error_on_warning" )
   endif()
endif ( )


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

if( HYTEG_DOWNLOAD_BOOST )
   include(FetchContent)
   FetchContent_Declare(
           boost
           URL      https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.bz2
           URL_HASH SHA256=af57be25cb4c4f4b413ed692fe378affb4352ea50fbe294a11ef548f4d527d89
   )
   FetchContent_MakeAvailable( boost )
   set(BOOST_ROOT ${hyteg_BINARY_DIR}/_deps/boost-src/ CACHE PATH "Path to boost root")
endif()

if( HYTEG_BUILD_WITH_PYTHON3 )
   find_package (Python3 COMPONENTS Interpreter Development)
endif()

# Add flags and libraries required to build terraneo submodule
if( HYTEG_TERRANEO_MODULE )


   find_package( Boost 1.48 )
   if ( NOT Boost_FOUND )
      message( FATAL_ERROR "Boost is required for geometric calculations on the sphere but was not found.
      Please define BOOST_ROOT (e.g. cmake -DBOOST_ROOT=/local/boost/)
      or enable the cmake option HYTEG_DOWNLOAD_BOOST (e.g. cmake -DHYTEG_DOWNLOAD_BOOST=ON)" )
   endif ()


   # Unpack plate topology file
   if( NOT EXISTS ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson )
     if( EXISTS ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson.bz2 )
       message( STATUS "Unpacking plate topologies file" )
       execute_process(COMMAND bunzip2 "--keep" "topologies0-100Ma.geojson.bz2" WORKING_DIRECTORY ${hyteg_SOURCE_DIR}/data/terraneo/plates )
     else()
       message( WARNING "Missing plate topologies file: ${hyteg_SOURCE_DIR}/data/terraneo/plates/topologies0-100Ma.geojson.bzip2" )
     endif()
   endif()

   # For using signed int in ADIOS2 output for post-processing with Fortran codes
   if( HYTEG_BUILD_WITH_ADIOS2 AND NOT HYTEG_USE_SIGNED_INT_FOR_ADIOS2 )
     set( HYTEG_USE_SIGNED_INT_FOR_ADIOS2 ON CACHE BOOL "Use int64_t when writing BP files for ParaView" FORCE )
     message( WARNING "HYTEG_USE_SIGNED_INT_FOR_ADIOS2 was forced ON as HYTEG_TERRANEO_MODULE=ON" )
   endif()

endif()

# Add flags and libraries required to build the mantle convection app
if( HYTEG_MANTLECONVECTION_APP )
   find_package( Boost 1.48 )
   if ( NOT Boost_FOUND )
      message( FATAL_ERROR "Boost is required for geometric calculations on the sphere but was not found.
      Please define BOOST_ROOT (e.g. cmake -DBOOST_ROOT=/local/boost/)
      or enable the cmake option HYTEG_DOWNLOAD_BOOST (e.g. cmake -DHYTEG_DOWNLOAD_BOOST=ON)" )
   endif ()
endif()

if ( HYTEG_BUILD_WITH_PETSC )
    find_package(PETSc REQUIRED)
    add_library(PETSc::PETSc INTERFACE IMPORTED)
    set_target_properties(PETSc::PETSc PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDES}"
          INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES}"
          )
    set(HYTEG_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE)
    set( WALBERLA_BUILD_WITH_MPI ON CACHE BOOL "Build with MPI" FORCE )
    message(STATUS "HYTEG_BUILD_WITH_MPI was force set to ON")
    set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDES})
    include(CheckSymbolExists)
    unset(HYTEG_PETSC_BUILT_WITH_HDF5 CACHE)
    check_symbol_exists(PETSC_HAVE_HDF5 petscconf.h  HYTEG_PETSC_BUILT_WITH_HDF5)
    if ( HYTEG_PETSC_BUILT_WITH_HDF5 )
       find_path( HDF5_PATH hdf5.h )
       if ( HDF5_PATH )
          message( STATUS ${HDF5_PATH} )
          target_include_directories( PETSc::PETSc SYSTEM INTERFACE "${HDF5_PATH}" )
       else ()
          message( WARNING "PETSc was build with hdf5 support but hdf5.h could not be found.\
          Add the path to hdf5.h with e.g. -DCMAKE_PREFIX_PATH=\${HDF5_INCLUDE_DIR}" )
       endif ()
    endif ()
endif()

if ( HYTEG_BUILD_WITH_TRILINOS )
   find_package( Trilinos REQUIRED )
   message(STATUS "Found Trilinos! Trilinos_DIR = ${Trilinos_DIR} ")
   INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
   list ( APPEND SERVICE_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARY_DIRS})
endif()

if ( HYTEG_BUILD_WITH_MPI )
   find_package(MPI REQUIRED)
endif ()

if ( HYTEG_BUILD_WITH_OPENMP )
   find_package(OpenMP REQUIRED)
endif ()

if ( HYTEG_BUILD_WITH_ADIOS2 )
  if ( HYTEG_BUILD_WITH_MPI )
    find_package( ADIOS2 REQUIRED COMPONENTS CXX MPI )
  else()
    find_package( ADIOS2 REQUIRED COMPONENTS CXX )
  endif()

  if ( ${ADIOS2_VERSION} VERSION_LESS "2.9.2" )
    message( WARNING "ADIOS2 version ${ADIOS2_VERSION} is below 2.9.2 !\nThus, your version is susceptible to the hashing issue described in\n https://github.com/ornladios/ADIOS2/discussions/3822\nThis can severly affect HyTeG's current checkpointing implementation!" )
  endif()
endif()

if ( HYTEG_BUILD_WITH_LIKWID )
   find_library( LIKWID_LIB likwid HINTS $ENV{LIKWID_LIBDIR} $ENV{LIKWID_ROOT}/lib )
   find_path( LIKWID_INCLUDE_DIR likwid.h HINTS $ENV{LIKWID_INCDIR} $ENV{LIKWID_ROOT}/include )


   if ( LIKWID_LIB AND LIKWID_INCLUDE_DIR)
      add_definitions ( "-DLIKWID_PERFMON" )
   else()
      message(WARNING "likwid marker library not found. Set environment variable LIKWID_ROOT")
      set ( HYTEG_BUILD_WITH_LIKWID OFF CACHE BOOL "Compile in markers for likwid-perfctr" FORCE )
   endif()
endif()

configure_file ( ${hyteg_SOURCE_DIR}/src/hyteg/HytegDefinitions.in.hpp
                 src/hyteg/HytegDefinitions.hpp )

#include_directories( ${hyteg_BINARY_DIR}/src )

list(APPEND ALL_DIRS data/meshes/2D data/meshes/3D data/param output doc/images data/kerncraftMachineFiles data/terraneo
                     data/terraneo/plates data/terraneo/viscosityProfiles data/terraneo/viscosityProfiles/tests data/terraneo/radialProfiles )

foreach(DIR IN LISTS ALL_DIRS)
   file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${DIR}")
endforeach ()

# Determine if Gmsh executable is available for testing purposes
find_program( GMSH_EXECUTABLE gmsh )
if( GMSH_EXECUTABLE ) 
    message(STATUS "Found Gmsh executable at ${GMSH_EXECUTABLE}" )
endif()

waLBerla_link_files_to_builddir(data/meshes/2D/*.msh)
waLBerla_link_files_to_builddir(data/meshes/3D/*.msh)
waLBerla_link_files_to_builddir(data/param/*.prm)
waLBerla_link_files_to_builddir(data/kerncraftMachineFiles/*.yml)
waLBerla_link_files_to_builddir(doc/images/*.png)
waLBerla_link_files_to_builddir(data/terraneo/*.json)
waLBerla_link_files_to_builddir(data/terraneo/viscosityProfiles/*.json)
waLBerla_link_files_to_builddir(data/terraneo/viscosityProfiles/tests/*.json)
waLBerla_link_files_to_builddir(data/terraneo/plates/*.geojson)
waLBerla_link_files_to_builddir(data/terraneo/plates/*.rot)
waLBerla_link_files_to_builddir(data/terraneo/radialProfiles/*.json)

add_subdirectory( src )
add_subdirectory( apps )
add_subdirectory( tests )
add_subdirectory( tutorials )


############################################################################################################################
# Documentation Generation
#
# Build documentation using Doxygen (www.doxygen.org)
############################################################################################################################

  find_package ( Doxygen 1.10.0 OPTIONAL_COMPONENTS dot )
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
  list(APPEND ALL_OPTIONS BUILD_WITH_MPI BUILD_WITH_OPENMP BUILD_WITH_ADIOS2 BUILD_WITH_AVX BUILD_WITH_MPFR BUILD_WITH_PETSC BUILD_WITH_PYTHON3
                          BUILD_WITH_TRILINOS USE_GENERATED_KERNELS TERRANEO_MODULE )


  foreach(OPTION IN LISTS ALL_OPTIONS)
     if( HYTEG_${OPTION} )
        set( DOXYGEN_HYTEG_${OPTION} "HYTEG_${OPTION}" )
     endif()
  endforeach ( )

  # Avoid problems on systems that have doxygen installed, but not graphviz/dot;
  # configure_file below will set Doxygen's HAVE_DOT to our SYSTEM_HAS_DOT
  set ( SYSTEM_HAS_DOT "NO" )
  if ( TARGET Doxygen::dot )
    set ( SYSTEM_HAS_DOT "YES")
  endif()

  configure_file ( ${hyteg_SOURCE_DIR}/doc/doxygen.config ${hyteg_BINARY_DIR}/doc/doxygen.cfg @ONLY )

  add_custom_target ( doc_hyteg   ${DOXYGEN_EXECUTABLE} ${hyteg_BINARY_DIR}/doc/doxygen.cfg
    COMMENT "Generating API documentation with Doxygen" VERBATIM )

  if ( ${DOXYGEN_VERSION} VERSION_LESS "1.11.0" )
    message( WARNING "Doxygen version is older than 1.11.0 !\nThus, your version does not include a patch for "
      "issue https://github.com/doxygen/doxygen/issues/10517\nwhich affects some of HyTeG's tutorials" )
  endif()

else ()

  message( STATUS "No (suitable version of) Doxygen was found. You might not be able to build the documentation (properly)." )

endif ( )


############################################################################################################################


############################################################################################################################
# Git commit hash and submodules
############################################################################################################################

find_package(Git QUIET REQUIRED)
if(EXISTS "${PROJECT_SOURCE_DIR}/.git")

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

   include(cmake/CheckGit.cmake)
   CheckGitSetup()

else()
   message(FATAL_ERROR "Project folder seems to be no Git repository? Missing ${PROJECT_SOURCE_DIR}/.git")
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
        ${CMAKE_CURRENT_SOURCE_DIR}/src/hyteg/BuildInfo.in.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/src/hyteg/BuildInfo.hpp
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
