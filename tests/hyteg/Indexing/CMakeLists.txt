waLBerla_add_test_executable( CommonIndexingTest CommonIndexingTest.cpp )
target_link_libraries       ( CommonIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME CommonIndexingTest)

waLBerla_add_test_executable( IndexBasisConversionTest IndexBasisConversionTest.cpp )
target_link_libraries       ( IndexBasisConversionTest hyteg walberla::core )
waLBerla_execute_test(NAME IndexBasisConversionTest)

waLBerla_add_test_executable( P1EdgeIndexingTest P1EdgeIndexingTest.cpp )
target_link_libraries       ( P1EdgeIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_add_test_executable( P1FaceIndexingTest P1FaceIndexingTest.cpp )
target_link_libraries       ( P1FaceIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_add_test_executable( VertexDoFMacroCellIndexingTest VertexDoFMacroCellIndexingTest.cpp )
target_link_libraries       ( VertexDoFMacroCellIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME VertexDoFMacroCellIndexingTest)

waLBerla_add_test_executable( MicroCellIterationAndIndexingTest MicroCellIterationAndIndexingTest.cpp )
target_link_libraries       ( MicroCellIterationAndIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME MicroCellIterationAndIndexingTest)

waLBerla_add_test_executable( EdgeDoFOnMacroEdgeIndexingTest EdgeDoFOnMacroEdgeIndexingTest.cpp )
target_link_libraries       ( EdgeDoFOnMacroEdgeIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFOnMacroEdgeIndexingTest)

waLBerla_add_test_executable( EdgeDoFOnMacroFaceIndexingTest EdgeDoFOnMacroFaceIndexingTest.cpp )
target_link_libraries       ( EdgeDoFOnMacroFaceIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFOnMacroFaceIndexingTest)

waLBerla_add_test_executable( DGIndexingTest DGIndexingTest.cpp )
target_link_libraries       ( DGIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME DGIndexingTest)

waLBerla_add_test_executable( EdgeDoFOrientationConversionTest EdgeDoFOrientationConversionTest.cpp )
target_link_libraries       ( EdgeDoFOrientationConversionTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFOrientationConversionTest)

waLBerla_add_test_executable( MicroCellFromMicroVerticesTest MicroCellFromMicroVerticesTest.cpp )
target_link_libraries       ( MicroCellFromMicroVerticesTest hyteg walberla::core )
waLBerla_execute_test(NAME MicroCellFromMicroVerticesTest)

waLBerla_add_test_executable( MicroVertexIdxMacroRefinementTest MicroVertexIdxMacroRefinementTest.cpp )
target_link_libraries       ( MicroVertexIdxMacroRefinementTest hyteg walberla::core )
waLBerla_execute_test(NAME MicroVertexIdxMacroRefinementTest)

waLBerla_add_test_executable( VolumeDoFIndexingTest VolumeDoFIndexingTest.cpp )
target_link_libraries       ( VolumeDoFIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME VolumeDoFIndexingTest)

waLBerla_add_test_executable( MicroElementDoFTypificationTest MicroElementDoFTypificationTest.cpp )
target_link_libraries       ( MicroElementDoFTypificationTest hyteg walberla::core )
waLBerla_execute_test(NAME MicroElementDoFTypificationTest)

waLBerla_add_test_executable( CubeIndexingTest CubeIndexingTest.cpp )
target_link_libraries       ( CubeIndexingTest hyteg walberla::core )
waLBerla_execute_test(NAME CubeIndexingTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( CouplingCountTest CouplingCountTest.cpp )
    target_link_libraries       ( CouplingCountTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME CouplingCountTest COMMAND $<TARGET_FILE:CouplingCountTest>)
    waLBerla_execute_test(NAME ParallelCouplingCountTest3 COMMAND $<TARGET_FILE:CouplingCountTest> PROCESSES 3)
endif ()
