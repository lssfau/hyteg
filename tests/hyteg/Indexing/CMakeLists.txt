waLBerla_compile_test(FILES CommonIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME CommonIndexingTest)

waLBerla_compile_test(FILES IndexBasisConversionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IndexBasisConversionTest)

waLBerla_compile_test(FILES P1EdgeIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES P1FaceIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES VertexDoFMacroCellIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFMacroCellIndexingTest)

waLBerla_compile_test(FILES MicroCellIterationAndIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroCellIterationAndIndexingTest)

waLBerla_compile_test(FILES EdgeDoFOnMacroEdgeIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOnMacroEdgeIndexingTest)

waLBerla_compile_test(FILES EdgeDoFOnMacroFaceIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOnMacroFaceIndexingTest)

waLBerla_compile_test(FILES DGIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME DGIndexingTest)

waLBerla_compile_test(FILES EdgeDoFOrientationConversionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOrientationConversionTest)

waLBerla_compile_test(FILES MicroCellFromMicroVerticesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroCellFromMicroVerticesTest)

waLBerla_compile_test(FILES MicroVertexIdxMacroRefinementTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroVertexIdxMacroRefinementTest)

waLBerla_compile_test(FILES VolumeDoFIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VolumeDoFIndexingTest)

waLBerla_compile_test(FILES MicroElementDoFTypificationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroElementDoFTypificationTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES CouplingCountTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME CouplingCountTest COMMAND $<TARGET_FILE:CouplingCountTest>)
    waLBerla_execute_test(NAME ParallelCouplingCountTest3 COMMAND $<TARGET_FILE:CouplingCountTest> PROCESSES 3)
endif ()

