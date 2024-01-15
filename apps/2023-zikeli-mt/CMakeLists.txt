if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmark")
endif()

waLBerla_link_files_to_builddir("*.prm")

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
    waLBerla_add_executable( NAME test_Float16_support
            FILES test_Float16_support.cpp
            DEPENDS hyteg core)
    target_compile_features( test_Float16_support PUBLIC cxx_std_23)
endif ()

waLBerla_add_executable( NAME mpfr_print_version
        FILES mpfr_print_version.cpp
        DEPENDS MPFR::MPFR)
