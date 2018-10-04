waLBerla_link_files_to_builddir(*.msh)

## Convergence ##
waLBerla_compile_test(FILES convergence/P1CGConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CGConvergenceTest)

waLBerla_compile_test(FILES convergence/P1CG3DConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CG3DConvergenceTest)

waLBerla_compile_test(FILES convergence/P1GMGConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1GMGConvergenceTest)

waLBerla_compile_test(FILES convergence/P1GMG3DConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1GMG3DConvergenceTest)

waLBerla_compile_test(FILES convergence/P1StokesBlockCGTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1StokesBlockCGTest)

waLBerla_compile_test(FILES convergence/P1CGTensorVarcoeffConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CGTensorVarcoeffConvergenceTest)

waLBerla_compile_test(FILES convergence/P1CGVarcoeffConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CGVarcoeffConvergenceTest)

waLBerla_compile_test(FILES convergence/P1CGVarcoordsConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CGVarcoordsConvergenceTest)

waLBerla_compile_test(FILES convergence/P1GSConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1GSConvergenceTest)

waLBerla_compile_test(FILES convergence/P1MinResConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1MinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P1StokesMinResConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1StokesMinResConvergenceTest)

if(CMAKE_BUILD_TYPE MATCHES "Release")
  waLBerla_compile_test(FILES convergence/P1Stokes3DMinResConvergenceTest.cpp DEPENDS tinyhhg_core core)
  waLBerla_execute_test(NAME P1Stokes3DMinResConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P2CGConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2CGConvergenceTest)

waLBerla_compile_test(FILES convergence/P2SmoothConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2SmoothConvergenceTest)

waLBerla_compile_test(FILES convergence/P2JacobiConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2JacobiConvergenceTest)

waLBerla_compile_test(FILES convergence/P2MinResConvergenceTest.cpp DEPENDS tinyhhg_core core)
walberla_execute_test(NAME P2MinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P2GMGHRefinementConvergenceTest.cpp DEPENDS tinyhhg_core core)
walberla_execute_test(NAME P2GMGHRefinementConvergenceTest)

waLBerla_compile_test(FILES convergence/P2GMGPRefinementConvergenceTest.cpp DEPENDS tinyhhg_core core)
walberla_execute_test(NAME P2GMGPRefinementConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1StokesMinResConvergenceTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2P1StokesMinResConvergenceTest)

## Indexing ##

waLBerla_compile_test(FILES Indexing/CommonIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME CommonIndexingTest)

waLBerla_compile_test(FILES Indexing/IndexBasisConversionTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME IndexBasisConversionTest)

waLBerla_compile_test(FILES Indexing/P1EdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/P1FaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES Indexing/VertexDoFMacroCellIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFMacroCellIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroEdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFOnMacroEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroFaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFOnMacroFaceIndexingTest)

waLBerla_compile_test(FILES Indexing/DGIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME DGIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOrientationConversionTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFOrientationConversionTest)

if(CMAKE_BUILD_TYPE MATCHES "Release")
  if(NOT WALBERLA_CXX_COMPILER_IS_MSVC )
    waLBerla_compile_test(FILES Indexing/ConstexprIndexTest.cpp DEPENDS tinyhhg_core core)
    waLBerla_execute_test(NAME ConstexprIndexTest)
  endif()
endif()

## VertexDoFSpace ##

waLBerla_compile_test(FILES vertexdofspace/VertexDoFFunctionTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFFunctionTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFMacroCellPackInfoTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFMacroCellPackInfoTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFMemoryTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFMemoryTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFBasicFunctionsTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFBasicFunctionsTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFStencilAssemblyTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFStencilAssemblyTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFAdditiveCommunicationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFAdditiveCommunicationTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFLinearProlongation3DTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME VertexDoFLinearProlongation3DTest)

## EdgeDoFSpace ##

waLBerla_compile_test(FILES edgedofspace/EdgeDoFMemoryTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFMemoryTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFFunctionTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFFunctionTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFEnumerateTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFEnumerateTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFCommTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFCommTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFInterpolationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFInterpolationTest)

waLBerla_compile_test(FILES edgedofspace/EdgeDoFFunction3DTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFFunction3DTest)

## Operators ##

waLBerla_compile_test(FILES operators/EdgeDoFToVertexDoFOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFToVertexDoFOperatorTest)

waLBerla_compile_test(FILES operators/EdgeDoFToEdgeDoFOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME EdgeDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES operators/SimpleVertexDoFToEdgeDoFOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME SimpleVertexDoFToEdgeDoFOperatorTest)

#################### P1 ####################

waLBerla_compile_test(FILES P1/P1MassOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1MassOperatorTest )
waLBerla_execute_test(NAME P1MassOperatorTestMPI COMMAND $<TARGET_FILE:P1MassOperatorTest> PROCESSES 2 )

waLBerla_compile_test(FILES P1/P1ElementwiseMassOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1ElementwiseMassOperatorTest )

waLBerla_compile_test(FILES P1/P1ElementwisePolarMassOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1ElementwisePolarMassOperatorTest )

waLBerla_compile_test(FILES P1/P1SmoothingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1SmoothingTest )

waLBerla_compile_test(FILES P1/P1LaplaceOperator3DTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1LaplaceOperator3DTest)

waLBerla_compile_test(FILES P1/P1InterpolateTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES P1/P1CommTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES P1/P1DataHandlingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1DataHandlingTest)


############################################

## P2 ##

waLBerla_compile_test(FILES P2/P2GSTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2GSTest)

waLBerla_compile_test(FILES P2/P2SimpleSmoothTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2SimpleSmoothTest)

waLBerla_compile_test(FILES P2/P2P1TransferTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2P1TransferTest)

waLBerla_compile_test(FILES P2/P2VertexInterpolateTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2VertexInterpolateTest)

waLBerla_compile_test(FILES P2/P2RestrictTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2RestrictTest)

waLBerla_compile_test(FILES P2/P2ProlongateTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2ProlongateTest)

waLBerla_compile_test(FILES P2/P2BasicFunctionsTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2BasicFunctionsTest)

waLBerla_compile_test(FILES P2/P2FunctionTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2FunctionTest)

waLBerla_compile_test(FILES P2/P2ConstantMassOperatorTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2ConstantMassOperatorTest )
waLBerla_execute_test(NAME P2ConstantMassOperatorTestMPI COMMAND $<TARGET_FILE:P2ConstantMassOperatorTest> PROCESSES 2 )

waLBerla_compile_test(FILES P2/P2LaplaceOperator3DTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P2LaplaceOperator3DTest)

if( HHG_BUILD_WITH_PETSC )
  waLBerla_compile_test(FILES P2/P2PetscSolveTest.cpp DEPENDS tinyhhg_core core )
  waLBerla_execute_test(NAME P2PetscSolveTest1 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 1)
  waLBerla_execute_test(NAME P2PetscSolveTest2 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 2)
  waLBerla_execute_test(NAME P2PetscSolveTest3 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 8)
endif()

## convergence ##

waLBerla_compile_test(FILES convergence/P2GSConvergenceTest.cpp DEPENDS tinyhhg_core core )
waLBerla_execute_test(NAME P2GSConvergenceTest)

if( HHG_BUILD_WITH_PETSC )
  waLBerla_compile_test(FILES P2/P2PetscTest.cpp DEPENDS tinyhhg_core core)
  waLBerla_execute_test(NAME P2PetscTest1 COMMAND $<TARGET_FILE:P2PetscTest> )
  waLBerla_execute_test(NAME P2PetscTest2 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 2)
  waLBerla_execute_test(NAME P2PetscTest8 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 8)
endif()

## Other ##

waLBerla_compile_test(FILES DGInterpolateTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME DGInterpolateTest)

waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME vCycleTest1 COMMAND $<TARGET_FILE:vCycleTest> )

waLBerla_compile_test(FILES DGCommTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME DGCommTest)

waLBerla_compile_test(FILES DGEnumerateTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME DGEnumerateTest)

waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES LevelinfoTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME LevelinfoTest)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveStorageTest1 COMMAND $<TARGET_FILE:PrimitiveStorageTest> )
waLBerla_execute_test(NAME PrimitiveStorageTest3 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 3 )
waLBerla_execute_test(NAME PrimitiveStorageTest8 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 8 )

waLBerla_compile_test(FILES PrimitiveRelativeOrderingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveRelativeOrderingTest)

waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME BufferedCommunicationTest1 COMMAND $<TARGET_FILE:BufferedCommunicationTest> )
waLBerla_execute_test(NAME BufferedCommunicationTest3 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 3 )
waLBerla_execute_test(NAME BufferedCommunicationTest8 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 8 )

waLBerla_compile_test(FILES PrimitiveMigrationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveMigrationTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> )
waLBerla_execute_test(NAME PrimitiveMigrationTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 3 )
waLBerla_execute_test(NAME PrimitiveMigrationTest8 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 8 )

if( WALBERLA_BUILD_WITH_PARMETIS )
  waLBerla_compile_test(FILES ParMetisTestHHG.cpp DEPENDS tinyhhg_core core)
  waLBerla_execute_test(NAME ParMetisTestHHG1 COMMAND $<TARGET_FILE:ParMetisTestHHG> )
  waLBerla_execute_test(NAME ParMetisTestHHG3 COMMAND $<TARGET_FILE:ParMetisTestHHG> PROCESSES 3 )
  waLBerla_execute_test(NAME ParMetisTestHHG8 COMMAND $<TARGET_FILE:ParMetisTestHHG> PROCESSES 8 )
endif()

waLBerla_compile_test(FILES MeshInfoTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME MeshInfoTest1 COMMAND $<TARGET_FILE:MeshInfoTest> )
waLBerla_execute_test(NAME MeshInfoTest2 COMMAND $<TARGET_FILE:MeshInfoTest> PROCESSES 2 )
waLBerla_compile_test(FILES MeshGenTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME MeshInfoTest3 COMMAND $<TARGET_FILE:MeshGenTest> )
waLBerla_execute_test(NAME MeshInfoTest4 COMMAND $<TARGET_FILE:MeshGenTest> PROCESSES 2 )

waLBerla_compile_test(FILES FunctionSpaceDataTypesTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME FunctionSpaceDataTypesTest1 COMMAND $<TARGET_FILE:FunctionSpaceDataTypesTest> )

if( HHG_BUILD_WITH_PETSC )
  waLBerla_compile_test(FILES PetscTest.cpp DEPENDS tinyhhg_core core)
  waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest> )
  waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2 )
endif()

waLBerla_compile_test(FILES FunctionMemorySerializationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> )
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2 )
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8 )

# Blending
waLBerla_compile_test(FILES blending/P1StokesBlendingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1StokesBlendingTest1 COMMAND $<TARGET_FILE:P1StokesBlendingTest> )
waLBerla_execute_test(NAME P1StokesBlendingTest2 COMMAND $<TARGET_FILE:P1StokesBlendingTest> PROCESSES 2 )
waLBerla_execute_test(NAME P1StokesBlendingTest8 COMMAND $<TARGET_FILE:P1StokesBlendingTest> PROCESSES 8 )

if( HHG_BUILD_WITH_EIGEN )
  waLBerla_compile_test(FILES blending/P1StokesPolynomialBlendingTest.cpp DEPENDS tinyhhg_core core)
  waLBerla_execute_test(NAME P1StokesPolynomialBlendingTest1 COMMAND $<TARGET_FILE:P1StokesPolynomialBlendingTest> )
  waLBerla_execute_test(NAME P1StokesPolynomialBlendingTest2 COMMAND $<TARGET_FILE:P1StokesPolynomialBlendingTest> PROCESSES 2 )
  waLBerla_execute_test(NAME P1StokesPolynomialBlendingTest8 COMMAND $<TARGET_FILE:P1StokesPolynomialBlendingTest> PROCESSES 8 )
endif()

# Test some basic auxiliary stuff for application layer
waLBerla_compile_test(FILES FindMaxMinMagTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME FindMaxMinMagTest1 COMMAND $<TARGET_FILE:FindMaxMinMagTest> )
waLBerla_execute_test(NAME FindMaxMinMagTest2 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 2 )
