waLBerla_add_test_executable( VertexDoFMacroCellPackInfoTest VertexDoFMacroCellPackInfoTest.cpp )
target_link_libraries       ( VertexDoFMacroCellPackInfoTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFMacroCellPackInfoTest)

waLBerla_add_test_executable( VertexDoFMemoryTest VertexDoFMemoryTest.cpp )
target_link_libraries       ( VertexDoFMemoryTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFMemoryTest)

waLBerla_add_test_executable( VertexDoFBasicFunctionsTest VertexDoFBasicFunctionsTest.cpp )
target_link_libraries       ( VertexDoFBasicFunctionsTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFBasicFunctionsTest)

waLBerla_add_test_executable( VertexDoFStencilAssemblyTest VertexDoFStencilAssemblyTest.cpp )
target_link_libraries       ( VertexDoFStencilAssemblyTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFStencilAssemblyTest)

waLBerla_add_test_executable( VertexDoFAdditiveCommunicationTest VertexDoFAdditiveCommunicationTest.cpp )
target_link_libraries       ( VertexDoFAdditiveCommunicationTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFAdditiveCommunicationTest)

waLBerla_add_test_executable( VertexDoFLinearProlongation3DTest VertexDoFLinearProlongation3DTest.cpp )
target_link_libraries       ( VertexDoFLinearProlongation3DTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFLinearProlongation3DTest)

waLBerla_add_test_executable( VertexDoFQuadraticProlongation3DTest VertexDoFQuadraticProlongation3DTest.cpp )
target_link_libraries       ( VertexDoFQuadraticProlongation3DTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFQuadraticProlongation3DTest)

waLBerla_add_test_executable( VertexDoFFunctionTest VertexDoFFunctionTest.cpp )
target_link_libraries       ( VertexDoFFunctionTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME VertexDoFFunctionTest)

