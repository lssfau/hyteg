waLBerla_compile_test(FILES basicMixedPrecisionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME basicMixedPrecisionTest)

waLBerla_compile_test(FILES mixedPrecisionIterativeRefinement.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME mixedPrecisionIterativeRefinement)

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
    waLBerla_compile_test(FILES float16SupportTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME float16SupportTest)
endif ()

