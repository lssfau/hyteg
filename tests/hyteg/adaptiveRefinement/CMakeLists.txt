waLBerla_add_test_executable( adaptiveRefinementTest adaptiveRefinementTest.cpp )
target_link_libraries       ( adaptiveRefinementTest hyteg walberla::core )
waLBerla_execute_test(NAME adaptiveRefinementTest)

waLBerla_add_test_executable( adaptiveRefinementBoundaryTest adaptiveRefinementBoundaryTest.cpp )
target_link_libraries       ( adaptiveRefinementBoundaryTest hyteg walberla::core )
waLBerla_execute_test(NAME adaptiveRefinementBoundaryTest)

if(WALBERLA_DOUBLE_ACCURACY)
waLBerla_add_test_executable( adaptiveRefinementBlendingTest adaptiveRefinementBlendingTest.cpp )
target_link_libraries       ( adaptiveRefinementBlendingTest hyteg walberla::core )
waLBerla_execute_test(NAME adaptiveRefinementBlendingTest)
endif()

