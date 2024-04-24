if ( HYTEG_BUILD_WITH_MPFR )
    waLBerla_compile_test(FILES buildWithMPFRTest.cpp DEPENDS core)
    waLBerla_execute_test(NAME buildWithMPFRTest)

    waLBerla_compile_test(FILES mpfr/MPFRWrapperTest.cpp DEPENDS core)
    waLBerla_execute_test(NAME MPFRWrapperTest)
endif ()
