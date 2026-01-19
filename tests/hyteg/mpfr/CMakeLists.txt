if ( HYTEG_BUILD_WITH_MPFR )
    waLBerla_add_test_executable( buildWithMPFRTest buildWithMPFRTest.cpp )
    target_link_libraries       ( buildWithMPFRTest MPFR::MPFR walberla::core )
    waLBerla_execute_test(NAME buildWithMPFRTest)

    waLBerla_add_test_executable( MPFRWrapperTest MPFRWrapperTest.cpp )
    target_link_libraries       ( MPFRWrapperTest MPFR::MPFR walberla::core )
    waLBerla_execute_test(NAME MPFRWrapperTest)
endif ()
