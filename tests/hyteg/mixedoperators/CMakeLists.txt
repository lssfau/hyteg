if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES P1ToDGOperatorTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1ToDGOperatorTest)
endif()

waLBerla_compile_test(FILES MixedDG1AndP2PlusBubbleOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MixedDG1AndP2PlusBubbleOperatorTest)
