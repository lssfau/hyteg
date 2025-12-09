if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmark")
endif()

waLBerla_link_files_to_builddir("*.prm")

if( HYTEG_BUILD_WITH_MPFR )
    add_executable       ( mpfr_print_version mpfr_print_version.cpp )
    target_link_libraries( mpfr_print_version MPFR::MPFR )
endif()

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
    add_subdirectory(MT-apps)
endif ()
