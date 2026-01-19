waLBerla_add_test_executable( MixedDG1AndP2PlusBubbleFormsTest MixedDG1AndP2PlusBubbleFormsTest.cpp )
target_link_libraries       ( MixedDG1AndP2PlusBubbleFormsTest hyteg walberla::core )
waLBerla_execute_test(NAME MixedDG1AndP2PlusBubbleFormsTest)

waLBerla_add_test_executable( P2LinearCombinationFormTest P2LinearCombinationFormTest.cpp )
target_link_libraries       ( P2LinearCombinationFormTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2LinearCombinationFormTest)

waLBerla_add_test_executable( RowSumTest RowSumTest.cpp )
target_link_libraries       ( RowSumTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME RowSumTest)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_add_test_executable( HytegVsFenicsFormTest HytegVsFenicsFormTest.cpp )
    target_link_libraries       ( HytegVsFenicsFormTest hyteg walberla::core )
    waLBerla_execute_test(NAME HytegVsFenicsFormTest)
endif()

waLBerla_add_test_executable( SingleRowIntegrationTest SingleRowIntegrationTest.cpp )
target_link_libraries       ( SingleRowIntegrationTest hyteg walberla::core )
waLBerla_execute_test(NAME SingleRowIntegrationTest)

