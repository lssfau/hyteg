if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmark")
endif()

waLBerla_link_files_to_builddir("*.prm")

if( HYTEG_BUILD_WITH_MPFR )
    waLBerla_add_executable( NAME mpfr_print_version
            FILES mpfr_print_version.cpp
            DEPENDS MPFR::MPFR)
endif()

add_subdirectory(MT-apps)
