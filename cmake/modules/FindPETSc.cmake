# - Try to find PETSc
# Once done this will define
#  PETSC_FOUND        - System has PETSc
#  PETSC_INCLUDE_DIRS - The PETSc include directories
#  PETSC_LIBRARY_DIR  - The PETSc library directory
#  PETSC_LIBRARIES    - The libraries needed to use PETSc
#
# Setting these changes the behavior of the search
#  PETSC_DIR  - directory in which PETSc resides
#  PETSC_ARCH - build architecture

set(PETSC_DIR "" CACHE PATH "An optional hint to a PETSc directory")
set(PETSC_ARCH "" CACHE STRING "An optional hint to a PETSc arch")

if("${PETSC_DIR}" STREQUAL "")
    set(PETSC_DIR "$ENV{PETSC_DIR}")
endif()


if("${PETSC_ARCH}" STREQUAL "")
    set(PETSC_ARCH "$ENV{PETSC_ARCH}")
endif()

set( ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig" )


find_package(PkgConfig)
pkg_check_modules(PC_PETSC PETSc)

set(PETSC_INCLUDE_DIRS ${PC_PETSC_INCLUDEDIR})
set(PETSC_LIBRARY_DIR ${PC_PETSC_STATIC_LIBRARY_DIRS})
set(PETSC_LIBRARIES ${PC_PETSC_STATIC_LIBRARIES})
set(PETSC_VERSION ${PC_PETSC_VERSION})


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PETSc "Could NOT find PETSc. Please make sure to set PETSC_DIR and PETSC_ARCH variables"
        PETSC_VERSION PETSC_LIBRARY_DIR PETSC_LIBRARIES PETSC_INCLUDE_DIRS)