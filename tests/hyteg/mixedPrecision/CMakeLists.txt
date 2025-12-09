waLBerla_add_test_executable( basicMixedPrecisionTest basicMixedPrecisionTest.cpp )
target_link_libraries       ( basicMixedPrecisionTest hyteg walberla::core )
waLBerla_execute_test(NAME basicMixedPrecisionTest)

waLBerla_add_test_executable( mixedPrecisionIterativeRefinement mixedPrecisionIterativeRefinement.cpp )
target_link_libraries       ( mixedPrecisionIterativeRefinement hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME mixedPrecisionIterativeRefinement)

if ( WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT )
    waLBerla_add_test_executable( float16SupportTest float16SupportTest.cpp )
    target_link_libraries       ( float16SupportTest hyteg walberla::core )
    waLBerla_execute_test(NAME float16SupportTest)
endif ()

