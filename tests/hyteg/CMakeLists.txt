waLBerla_link_files_to_builddir(*.msh)

## adaptive refinement ##
waLBerla_compile_test(FILES adaptiveRefinement/adaptiveRefinementTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME adaptiveRefinementTest)

waLBerla_compile_test(FILES adaptiveRefinement/adaptiveRefinementBoundaryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME adaptiveRefinementBoundaryTest)

if(WALBERLA_DOUBLE_ACCURACY)
waLBerla_compile_test(FILES adaptiveRefinement/adaptiveRefinementBlendingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME adaptiveRefinementBlendingTest)
endif()

waLBerla_compile_test(FILES convergence/ElementwiseDivKGradCGConvergenceTest.cpp DEPENDS hyteg core opgen-div_k_grad)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradCGConvergenceTest)
endif ()

waLBerla_compile_test(FILES convergence/ElementwiseDivKGradNeumannCGConvergenceTest.cpp DEPENDS hyteg core opgen-div_k_grad)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradNeumannCGConvergenceTest)
endif ()

## convection ##

## dataexport+import ##

if (HYTEG_BUILD_WITH_ADIOS2)
  waLBerla_compile_test(FILES dataexport+import/BasicADIOSTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME BasicADIOSTest)

  waLBerla_compile_test(FILES dataexport+import/AdiosWriterTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME AdiosWriterTest)

  waLBerla_compile_test(FILES dataexport+import/CheckpointingTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME CheckpointingTest)

  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/dataexport+import/CheckpointRestoreTest.prm
    ${CMAKE_CURRENT_BINARY_DIR}/CheckpointRestoreTest.prm )
  waLBerla_compile_test(FILES dataexport+import/CheckpointRestoreTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME CheckpointRestoreTest-Sequential COMMAND $<TARGET_FILE:CheckpointRestoreTest> PROCESSES 1)
  waLBerla_execute_test(NAME CheckpointRestoreTest-Parallel-Same-Communicator COMMAND $<TARGET_FILE:CheckpointRestoreTest> PROCESSES 3)
  waLBerla_execute_test(NAME CheckpointRestoreTest-Parallel-Import-Only COMMAND $<TARGET_FILE:CheckpointRestoreTest> CheckpointRestoreTest.prm -Parameters.onlyImport=true PROCESSES 3)
endif ()

waLBerla_compile_test(FILES dataexport+import/TableTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME TableTest)
set_tests_properties(TableTest PROPERTIES FIXTURES_SETUP TableTest)

if(WALBERLA_DOUBLE_ACCURACY)
    set(TABLE_REFERENCE_FILE "TableTest.dat.ref.dp")
else()
    set(TABLE_REFERENCE_FILE "TableTest.dat.ref.sp")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/dataexport+import/${TABLE_REFERENCE_FILE}
      ${CMAKE_CURRENT_BINARY_DIR}/${TABLE_REFERENCE_FILE})
ADD_TEST(TableTestCheckOutput ${CMAKE_COMMAND} -E compare_files TableTest.dat ${TABLE_REFERENCE_FILE})
set_tests_properties(TableTestCheckOutput PROPERTIES FIXTURES_REQUIRED TableTest)

waLBerla_compile_test(FILES dataexport+import/KeyValueStoreTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME KeyValueStoreTest)
set_tests_properties(KeyValueStoreTest PROPERTIES FIXTURES_SETUP KeyValueStoreTest)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/dataexport+import/KeyValueStoreTest.tex.ref
      ${CMAKE_CURRENT_BINARY_DIR}/KeyValueStoreTest.tex.ref)
ADD_TEST(KeyValueStoreTestCheckOutput ${CMAKE_COMMAND} -E compare_files KeyValueStoreTest.tex KeyValueStoreTest.tex.ref)
set_tests_properties(KeyValueStoreTestCheckOutput PROPERTIES FIXTURES_REQUIRED KeyValueStoreTest)

## Forms ##

## Indexing ##

waLBerla_compile_test(FILES Indexing/CommonIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME CommonIndexingTest)

waLBerla_compile_test(FILES Indexing/IndexBasisConversionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IndexBasisConversionTest)

waLBerla_compile_test(FILES Indexing/P1EdgeIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/P1FaceIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES Indexing/VertexDoFMacroCellIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFMacroCellIndexingTest)

waLBerla_compile_test(FILES Indexing/MicroCellIterationAndIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroCellIterationAndIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroEdgeIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOnMacroEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroFaceIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOnMacroFaceIndexingTest)

waLBerla_compile_test(FILES Indexing/DGIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME DGIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOrientationConversionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFOrientationConversionTest)

waLBerla_compile_test(FILES Indexing/MicroCellFromMicroVerticesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroCellFromMicroVerticesTest)

waLBerla_compile_test(FILES Indexing/MicroVertexIdxMacroRefinementTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MicroVertexIdxMacroRefinementTest)

waLBerla_compile_test(FILES Indexing/VolumeDoFIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VolumeDoFIndexingTest)

## Math ##

waLBerla_compile_test(FILES math/IntersectionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IntersectionTest)

waLBerla_compile_test(FILES meshes/RefineCoarseMeshComparisonTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME RefineCoarseMeshComparisonTest)

## VertexDoFSpace ##

waLBerla_compile_test(FILES vertexdofspace/VertexDoFMacroCellPackInfoTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFMacroCellPackInfoTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFMemoryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFMemoryTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFBasicFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFBasicFunctionsTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFStencilAssemblyTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFStencilAssemblyTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFAdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFAdditiveCommunicationTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFLinearProlongation3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFLinearProlongation3DTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFQuadraticProlongation3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFQuadraticProlongation3DTest)

## EdgeDoFSpace ##

waLBerla_compile_test(FILES edgedofspace/EdgeDoFAdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFAdditiveCommunicationTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFMemoryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFMemoryTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFFunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFFunctionTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFEnumerateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFEnumerateTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFCommTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFCommTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFInterpolationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFInterpolationTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFFunction3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFFunction3DTest)

## Operators ##

waLBerla_compile_test(FILES operators/EdgeDoFToVertexDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFToVertexDoFOperatorTest)

waLBerla_compile_test(FILES operators/SimpleVertexDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SimpleVertexDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES operators/P1SurrogateOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1SurrogateOperatorTest)

#################### P0 ####################

waLBerla_compile_test(FILES P0/P0FunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P0FunctionTest1 COMMAND $<TARGET_FILE:P0FunctionTest>)
waLBerla_execute_test(NAME P0FunctionTest2 COMMAND $<TARGET_FILE:P0FunctionTest> PROCESSES 2)

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES mixedoperators/P1ToDGOperatorTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1ToDGOperatorTest)
endif()

#if( HYTEG_BUILD_WITH_PETSC )
#    waLBerla_compile_test(FILES dg/P1WithDGFromPoisson2DConvergenceTest.cpp DEPENDS hyteg core)
#    waLBerla_execute_test(NAME P1WithDGFromPoisson2DConvergenceTest)
#endif()

#################### P1 ####################

waLBerla_compile_test(FILES P1/P1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES P1/P1CommTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES P1/P1DataHandlingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1DataHandlingTest)

waLBerla_compile_test(FILES P1/P1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1EvaluateTest)

############################################

## P2 ##

waLBerla_compile_test(FILES P2/P2P1TransferTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1TransferTest)

waLBerla_compile_test(FILES P2/P2VertexInterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2VertexInterpolateTest)

waLBerla_compile_test(FILES P2/P2RestrictTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2RestrictTest)

waLBerla_compile_test(FILES P2/P2ProlongateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2ProlongateTest)

waLBerla_compile_test(FILES P2/P2QuadraticRestriction3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2QuadraticRestriction3DTest)

waLBerla_compile_test(FILES P2/P2QuadraticProlongation3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2QuadraticProlongation3DTest)

waLBerla_compile_test(FILES P2/P2BasicFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2BasicFunctionsTest)

waLBerla_compile_test(FILES P2/P2SwapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2SwapTest)

waLBerla_compile_test(FILES P2/P2FunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2FunctionTest)

waLBerla_compile_test(FILES P2/P2Level0Test.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2Level0Test1 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 1)
waLBerla_execute_test(NAME P2Level0Test2 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 2)

waLBerla_compile_test(FILES P2/P2EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2EvaluateTest)

## composites ##

waLBerla_compile_test(FILES composites/P0P1HelperFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest1 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 1)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest3 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 3)

if (HYTEG_BUILD_WITH_PETSC)
    #waLBerla_compile_test(FILES composites/P1P0Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME P1P0Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P0Stokes2DPetscSolveTest> PROCESSES 1)
endif ()

############################################

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES VectorFunctions/VectorFunctionPetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME VectorFunctionPetscTest)
endif ()

## Communication ##
waLBerla_compile_test(FILES communication/BufferedCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BufferedCommunicationTest1 COMMAND $<TARGET_FILE:BufferedCommunicationTest>)
waLBerla_execute_test(NAME BufferedCommunicationTest3 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 3)
waLBerla_execute_test(NAME BufferedCommunicationTest8 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 8)

waLBerla_compile_test(FILES communication/SyncingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SyncingTest1 COMMAND $<TARGET_FILE:SyncingTest>)
waLBerla_execute_test(NAME SyncingTest3 COMMAND $<TARGET_FILE:SyncingTest> PROCESSES 3)

## Other ##




waLBerla_compile_test(FILES primitives/PrimitiveID+TypeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveID+TypeTest)

waLBerla_compile_test(FILES primitives/PrimitiveDataTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES LevelinfoTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME LevelinfoTest)

waLBerla_compile_test(FILES primitives/rootOnlySetupStorageTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME rootOnlySetupStorageTest COMMAND $<TARGET_FILE:rootOnlySetupStorageTest>)
waLBerla_execute_test(NAME rootOnlySetupStorageTest3 COMMAND $<TARGET_FILE:rootOnlySetupStorageTest> PROCESSES 3)
waLBerla_execute_test(NAME rootOnlySetupStorageTest8 COMMAND $<TARGET_FILE:rootOnlySetupStorageTest> PROCESSES 8)

waLBerla_compile_test(FILES primitives/PrimitiveStorageTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageTest1 COMMAND $<TARGET_FILE:PrimitiveStorageTest>)
waLBerla_execute_test(NAME PrimitiveStorageTest3 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageTest8 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 8)

waLBerla_compile_test(FILES primitives/PrimitiveRelativeOrderingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveRelativeOrderingTest)

if (WALBERLA_BUILD_WITH_MPI)
    waLBerla_compile_test(FILES primitives/PrimitiveStorageParallelSetupTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest1 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest>)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest2 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 2)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest3 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 3)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest4 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 4)
endif()

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveMigrationTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationTest>)
waLBerla_execute_test(NAME PrimitiveMigrationTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationTest8 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/ParallelRoundRobinTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ParallelRoundRobinTest1 COMMAND $<TARGET_FILE:ParallelRoundRobinTest>)
waLBerla_execute_test(NAME ParallelRoundRobinTest3 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 3)
waLBerla_execute_test(NAME ParallelRoundRobinTest8 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/PrimitiveStorageCopyTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest1 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest>)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest3 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest8 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 8)

if (WALBERLA_BUILD_WITH_PARMETIS)
    waLBerla_compile_test(FILES adaptivity/ParMetisTestHyteg.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME ParMetisTestHyteg1 COMMAND $<TARGET_FILE:ParMetisTestHyteg>)
    waLBerla_execute_test(NAME ParMetisTestHyteg3 COMMAND $<TARGET_FILE:ParMetisTestHyteg> PROCESSES 3)
    waLBerla_execute_test(NAME ParMetisTestHyteg8 COMMAND $<TARGET_FILE:ParMetisTestHyteg> PROCESSES 8)
endif ()

waLBerla_compile_test(FILES PackageBufferSystemTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PackageBufferSystemTest1 COMMAND $<TARGET_FILE:PackageBufferSystemTest>)
waLBerla_execute_test(NAME PackageBufferSystemTest2 COMMAND $<TARGET_FILE:PackageBufferSystemTest> PROCESSES 2)
waLBerla_execute_test(NAME PackageBufferSystemTest3 COMMAND $<TARGET_FILE:PackageBufferSystemTest> PROCESSES 3)
waLBerla_execute_test(NAME PackageBufferSystemTest4 COMMAND $<TARGET_FILE:PackageBufferSystemTest> PROCESSES 4)

waLBerla_compile_test(FILES meshes/MeshInfoTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MeshInfoTest1 COMMAND $<TARGET_FILE:MeshInfoTest>)
waLBerla_execute_test(NAME MeshInfoTest2 COMMAND $<TARGET_FILE:MeshInfoTest> PROCESSES 2)
waLBerla_compile_test(FILES meshes/MeshGenTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MeshInfoTest3 COMMAND $<TARGET_FILE:MeshGenTest>)
waLBerla_execute_test(NAME MeshInfoTest4 COMMAND $<TARGET_FILE:MeshGenTest> PROCESSES 2)

waLBerla_compile_test(FILES meshes/2DCoarseMeshRefinementTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME 2DCoarseMeshRefinementTest)

waLBerla_compile_test(FILES functions/FunctionMemoryAllocationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest1 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest>)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest2 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest8 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 8)

waLBerla_compile_test(FILES functions/FunctionSpaceDataTypesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionSpaceDataTypesTest1 COMMAND $<TARGET_FILE:FunctionSpaceDataTypesTest>)

waLBerla_compile_test(FILES functions/FunctionIteratorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionIteratorTest1 COMMAND $<TARGET_FILE:FunctionIteratorTest>)
waLBerla_execute_test(NAME FunctionIteratorTest2 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionIteratorTest8 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 8)

waLBerla_compile_test(FILES functions/FunctionMultiStoreTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMultiStoreTest)

waLBerla_compile_test(FILES functions/FunctionRuleOfThreeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest1 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest>)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest2 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest8 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 8)

waLBerla_compile_test(FILES functions/FEFunctionRegistryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FEFunctionRegistryTest COMMAND $<TARGET_FILE:FEFunctionRegistryTest>)

## Blending ##
waLBerla_compile_test(FILES blending/AffineMap3DTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME AffineMap3DTest1 COMMAND $<TARGET_FILE:AffineMap3DTest>)
endif ()

waLBerla_compile_test(FILES blending/AnnulusMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME AnnulusMapTest1 COMMAND $<TARGET_FILE:AnnulusMapTest>)

waLBerla_compile_test(FILES blending/BasicBlendingMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BasicBlendingMapTest)

waLBerla_compile_test(FILES blending/IcosahedralShellMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IcosahedralShellMapTest COMMAND $<TARGET_FILE:IcosahedralShellMapTest>)

waLBerla_compile_test(FILES blending/InverseBlending+PrimitiveIDFindingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME InverseBlending+PrimitiveIDFindingTest)

waLBerla_compile_test(FILES blending/PolygonDistanceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PolygonDistanceTest)

waLBerla_compile_test(FILES blending/TokamakMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME TokamakMapTest)

# Test some basic auxiliary stuff for application layer
waLBerla_compile_test(FILES functions/FindMaxMinMagTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FindMaxMinMagTest1 COMMAND $<TARGET_FILE:FindMaxMinMagTest>)
waLBerla_execute_test(NAME FindMaxMinMagTest2 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 2)

waLBerla_compile_test(FILES functions/FunctionPropertiesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionPropertiesTest)

waLBerla_compile_test(FILES functions/FunctionMultElementwiseTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMultElementwiseTest)

waLBerla_compile_test(FILES functions/FunctionExtendedExpressionInterpolationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionExtendedExpressionInterpolationTest)

waLBerla_compile_test(FILES numerictools/L2Test.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME L2Test)
waLBerla_execute_test(NAME L2TestMPI COMMAND $<TARGET_FILE:L2Test> PROCESSES 3)

## Form Evaluation ##
if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES forms/HytegVsFenicsFormTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME HytegVsFenicsFormTest)
endif()

waLBerla_compile_test(FILES forms/SingleRowIntegrationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SingleRowIntegrationTest)

## Floating point exceptions ##
waLBerla_compile_test(FILES fpe/invertElementwiseTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME invertElementwiseTest)
waLBerla_execute_test(NAME invertElementwiseTestMPI COMMAND $<TARGET_FILE:invertElementwiseTest> PROCESSES 3)

waLBerla_compile_test(FILES fpe/CircularMapFPETest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME CircularMapFPETest)
waLBerla_execute_test(NAME CircularMapFPETestMPI COMMAND $<TARGET_FILE:CircularMapFPETest> PROCESSES 2)

waLBerla_compile_test(FILES freeslip/ProjectNormalTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ProjectNormalTest)
waLBerla_execute_test(NAME ProjectNormalTestMPI COMMAND $<TARGET_FILE:ProjectNormalTest> PROCESSES 2)

## Trilinos ##
if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES trilinos/TrilinosVersionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosVersionTest)
endif ()

waLBerla_compile_test(FILES meshes/TokamakLaplaceTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME TokamakLaplaceTest)
    waLBerla_execute_test(NAME TokamakLaplaceTestLongrun COMMAND $<TARGET_FILE:TokamakLaplaceTest> --longrun LABELS longrun)
endif ()

## BoundaryUID stuff ##
waLBerla_compile_test(FILES boundaryUID/BoundaryUIDTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BoundaryUIDTest)
waLBerla_compile_test(FILES boundaryUID/HollowBodyUIDTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME HollowBodyUIDTest)

## BlockFunctions + BlockOperators ##
waLBerla_compile_test(FILES BlockFunctions/FunctionWrapperTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionWrapperTest)
waLBerla_execute_test(NAME FunctionWrapperTestMPI COMMAND $<TARGET_FILE:FunctionWrapperTest> PROCESSES 2)

waLBerla_compile_test(FILES BlockFunctions/BlockFunctionBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BlockFunctionBasicTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Petsc/PetscCreateFromLinCombTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PetscCreateFromLinCombTest COMMAND $<TARGET_FILE:PetscCreateFromLinCombTest>)
endif ()

waLBerla_compile_test(FILES geometry/SphereToolsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SphereToolsTest1 COMMAND $<TARGET_FILE:SphereToolsTest>)
waLBerla_execute_test(NAME SphereToolsTest2 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 2)
waLBerla_execute_test(NAME SphereToolsTest4 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 4)

## DG ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES dg/DGPoisson2DConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 2 LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 4 LABELS longrun)

    waLBerla_compile_test(FILES dg/DGPoisson3DConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPoisson3DConvergenceTest)
    # Comment-in once macro-cell to macro-cell packing is implemented
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 2)
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 4)

    waLBerla_compile_test(FILES dg/DGLaplaceSymmetryTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGLaplaceSymmetryTest)

    waLBerla_compile_test(FILES dg/DGInterpolateEvaluateTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGInterpolateEvaluateTest)

    waLBerla_compile_test(FILES egfunctionspace/EGFunctionspaceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME EGFunctionspaceTest)


    waLBerla_compile_test(FILES dg/DGPetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPetscApplyTest)

    waLBerla_compile_test(FILES dg/DGGlobalRefinementGridConvTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGGlobalRefinementGridConvTest)

    waLBerla_compile_test(FILES dg/DGAdaptiveRefinementInterpolateEvaluateTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGAdaptiveRefinementInterpolateEvaluateTest)

    waLBerla_compile_test(FILES dg/DGAdaptiveRefinementGridConvTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGAdaptiveRefinementGridConvTest)

    waLBerla_compile_test(FILES dg/DGProlongationTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGProlongationTest)
endif ()

## Structural
waLBerla_compile_test(FILES PrintInfoTest.cpp)
waLBerla_execute_test(NAME PrintInfoTest)

## adaptive refinement w/ hanging nodes
waLBerla_compile_test(FILES adaptiverefinementhangingnodes/ARHangingNodesStorageTest.cpp)
waLBerla_execute_test(NAME ARHangingNodesStorageTest)


waLBerla_compile_test(FILES adaptiverefinementhangingnodes/MacroVolumeRefinementMicroIdxConversionTest.cpp)
waLBerla_execute_test(NAME MacroVolumeRefinementMicroIdxConversionTest)

waLBerla_compile_test(FILES eigen/EigenBasicTest.cpp)
waLBerla_execute_test(NAME EigenBasicTest)
waLBerla_execute_test(NAME EigenBasicTestMPI COMMAND $<TARGET_FILE:EigenBasicTest> PROCESSES 2)

waLBerla_compile_test(FILES eigen/EigenAssemblyTest.cpp)
waLBerla_execute_test(NAME EigenAssemblyTest)

waLBerla_compile_test(FILES eigen/EigenSparseSolverTest.cpp)
waLBerla_execute_test(NAME EigenSparseSolverTest)

## N1E1 ##
# For some reason executing some tests takes forever with clang and gcc in Debug mode.
set(N1E1_CONDITIONAL_LONGRUN "")
if (     CMAKE_BUILD_TYPE MATCHES "Debug"
     AND (WALBERLA_CXX_COMPILER_IS_CLANG OR WALBERLA_CXX_COMPILER_IS_GNU)
   )
    set(N1E1_CONDITIONAL_LONGRUN "longrun")
endif ()

waLBerla_compile_test(FILES N1E1/N1E1AddTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AddTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1AssignTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AssignTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1EvaluateTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1InterpolateTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1AdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AdditiveCommunicationTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1CurlCurlNullSpaceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlNullSpaceTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1MassTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1MassTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1FormCurlCurlTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1FormCurlCurlTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1CurlCurlConvergenceTest.cpp N1E1/common.hpp N1E1/common.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1CurlCurlConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

waLBerla_compile_test(FILES N1E1/N1E1ProlongationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1ProlongationTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1RestrictionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1RestrictionTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

waLBerla_compile_test(FILES N1E1/P1toN1E1GradientTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1toN1E1GradientTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1toP1LiftingTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1toP1LiftingTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1InverseDiagonalTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1InverseDiagonalTest LABELS n1e1)
endif ()

waLBerla_compile_test(FILES N1E1/N1E1BlendingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1BlendingTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1LinearFormBlendingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1LinearFormBlendingTest LABELS n1e1)

## types ##
waLBerla_compile_test(FILES types/MatrixConstructFromEigenTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MatrixConstructFromEigenTest)

waLBerla_compile_test(FILES mixedPrecision/basicMixedPrecisionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME basicMixedPrecisionTest)

## mpfr ##
if ( HYTEG_BUILD_WITH_MPFR )
    waLBerla_compile_test(FILES mpfr/buildWithMPFRTest.cpp DEPENDS core)
    waLBerla_execute_test(NAME buildWithMPFRTest)

    waLBerla_compile_test(FILES mpfr/MPFRWrapperTest.cpp DEPENDS core)
    waLBerla_execute_test(NAME MPFRWrapperTest)
endif ()
