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

## Convergence ##
waLBerla_compile_test(FILES convergence/P1CGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1CGConvergenceTest)
waLBerla_execute_test(NAME P1CGConvergenceTestMPI COMMAND $<TARGET_FILE:P1CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1CG3DConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1CG3DConvergenceTest)
waLBerla_execute_test(NAME P1CG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1CG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1GMGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1GMGConvergenceTest)
waLBerla_execute_test(NAME P1GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1FASConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1FASConvergenceTest)
waLBerla_execute_test(NAME P1FASConvergenceTestMPI COMMAND $<TARGET_FILE:P1FASConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1TauExtrapolationConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTest)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTestMPI COMMAND $<TARGET_FILE:P1TauExtrapolationConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1GMG3DConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1GMG3DConvergenceTest)
waLBerla_execute_test(NAME P1GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2GMGConvergenceTest)
waLBerla_execute_test(NAME P2GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMG3DConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2GMG3DConvergenceTest)
waLBerla_execute_test(NAME P2GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1VectorLaplaceCGTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1VectorLaplaceCGTest)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P1GSConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1GSConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P1JacobiConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1JacobiConvergenceTest)

waLBerla_compile_test(FILES convergence/P1MinResConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1MinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1StokesMinResConvergenceTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES convergence/P1P1StokesUzawaConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1P1StokesUzawaConvergenceTest)
endif ()

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_compile_test(FILES convergence/P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1Stokes3DMinResConvergenceTest)

    waLBerla_compile_test(FILES convergence/P2P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P1Stokes3DMinResConvergenceTest)

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES convergence/P1P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core)
        waLBerla_execute_test(NAME P1P1Stokes3DUzawaConvergenceTest)
    endif ()

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES convergence/P2P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core)
        waLBerla_execute_test(NAME P2P1Stokes3DUzawaConvergenceTest)
    endif ()
endif ()

waLBerla_compile_test(FILES convergence/P2CGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2CGConvergenceTest)
waLBerla_execute_test(NAME P2CGConvergenceTestMPI COMMAND $<TARGET_FILE:P2CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P1ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P2ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/ElementwiseDivKGradCGConvergenceTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradCGConvergenceTest)
endif ()

waLBerla_compile_test(FILES convergence/ElementwiseEpsilonMinResConvergenceTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTest)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTestLongrun COMMAND $<TARGET_FILE:ElementwiseEpsilonMinResConvergenceTest> --longrun LABELS longrun)
    endif ()
endif ()

waLBerla_compile_test(FILES convergence/P2SmoothConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2SmoothConvergenceTest)

waLBerla_compile_test(FILES convergence/P23DSmoothConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P23DSmoothConvergenceTest)

waLBerla_compile_test(FILES convergence/P2JacobiConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2JacobiConvergenceTest)

waLBerla_compile_test(FILES convergence/P2MinResConvergenceTest.cpp DEPENDS hyteg core)
walberla_execute_test(NAME P2MinResConvergenceTest)
waLBerla_execute_test(NAME P2MinResConvergenceTestMPI COMMAND $<TARGET_FILE:P2MinResConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGHRefinementConvergenceTest.cpp DEPENDS hyteg core)
walberla_execute_test(NAME P2GMGHRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGHRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGHRefinementConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGPRefinementConvergenceTest.cpp DEPENDS hyteg core)
walberla_execute_test(NAME P2GMGPRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGPRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGPRefinementConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1StokesMinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1StokesSchurCGConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1StokesSchurCGConvergenceTest)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P2GSConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2GSConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P2P1ElementwiseUzawaConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1ElementwiseUzawaConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1UzawaConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1UzawaConvergenceTest)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P1ChebyshevSmootherConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1ChebyshevSmootherConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P2UnsteadyDiffusion2DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2UnsteadyDiffusion2DTest)

waLBerla_compile_test(FILES convergence/P2UnsteadyDiffusionBlending2DTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME P2UnsteadyDiffusionBlending2DTest)
    endif ()
endif ()

waLBerla_compile_test(FILES convergence/StokesVarViscPopovTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME StokesVarViscPopovTest)
    endif ()
endif ()

## convection ##

waLBerla_compile_test(FILES convection/Upwind2DCircularConvectionTest.cpp DEPENDS hyteg core)
waLBerla_compile_test(FILES convection/AlgebraicUpwind3DCircularConvectionTest.cpp DEPENDS hyteg core)

## dataexport ##

waLBerla_compile_test(FILES dataexport/VTKOutputTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VTKOutputTest)

## Forms ##

waLBerla_compile_test(FILES forms/P2LinearCombinationFormTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2LinearCombinationFormTest)

waLBerla_compile_test(FILES forms/RowSumTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME RowSumTest)

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

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Indexing/CouplingCountTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME CouplingCountTest COMMAND $<TARGET_FILE:CouplingCountTest>)
    waLBerla_execute_test(NAME ParallelCouplingCountTest3 COMMAND $<TARGET_FILE:CouplingCountTest> PROCESSES 3)
endif ()

waLBerla_compile_test(FILES Indexing/VolumeDoFIndexingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VolumeDoFIndexingTest)

## Math ##

waLBerla_compile_test(FILES math/IntersectionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IntersectionTest)

## meshes ##
waLBerla_compile_test(FILES meshes/RefineCoarseMeshCGTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME RefineCoarseMeshCGTest)

waLBerla_compile_test(FILES meshes/RefineCoarseMeshComparisonTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME RefineCoarseMeshComparisonTest)

## VertexDoFSpace ##

waLBerla_compile_test(FILES vertexdofspace/VertexDoFFunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VertexDoFFunctionTest)

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

waLBerla_compile_test(FILES operators/EdgeDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES operators/SimpleVertexDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SimpleVertexDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES operators/MassOperatorTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME MassOperatorTest)
    waLBerla_execute_test(NAME MassOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 2)
endif ()

waLBerla_compile_test(FILES operators/DivergenceOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME DivergenceOperatorTest)

waLBerla_compile_test(FILES operators/ElementwiseOperatorAdditiveApplyTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ElementwiseOperatorAdditiveApplyTest)

waLBerla_compile_test(FILES operators/DiagonalNonConstantOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME DiagonalNonConstantOperatorTest)

waLBerla_compile_test(FILES operators/P1SurrogateOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1SurrogateOperatorTest)

waLBerla_compile_test(FILES operators/ViscousOperatorsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ViscousOperatorsTest)

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

waLBerla_compile_test(FILES P1/P1LaplaceOperator3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1LaplaceOperator3DTest)

waLBerla_compile_test(FILES P1/P1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES P1/P1CommTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES P1/P1DataHandlingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1DataHandlingTest)

waLBerla_compile_test(FILES P1/P1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1EvaluateTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1/P1PetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1PetscApplyTest1 COMMAND $<TARGET_FILE:P1PetscApplyTest>)
    waLBerla_execute_test(NAME P1PetscApplyTest2 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscApplyTest8 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1/P1PetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1PetscSolveTest1 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1PetscSolveTest2 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscSolveTest3 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 8)
endif ()

############################################

## P2 ##

waLBerla_compile_test(FILES P2/P2GSTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2GSTest)

waLBerla_compile_test(FILES P2/P2SimpleSmoothTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2SimpleSmoothTest)

waLBerla_compile_test(FILES P2/P2JacobiGSSmoothTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2JacobiGSSmoothTest)

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

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2LaplacePetscSymmetryTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2LaplacePetscSymmetryTest)
endif ()

waLBerla_compile_test(FILES P2/P2LaplaceOperator3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2PetscSolveTest1 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2PetscSolveTest2 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 2)
    if (NOT WALBERLA_CXX_COMPILER_IS_INTELLLVM)
        waLBerla_execute_test(NAME P2PetscSolveTest3 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 8)
    endif ()
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2PetscTest1 COMMAND $<TARGET_FILE:P2PetscTest>)
    waLBerla_execute_test(NAME P2PetscTest2 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscTest8 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2PetscApplyTest1 COMMAND $<TARGET_FILE:P2PetscApplyTest>)
    waLBerla_execute_test(NAME P2PetscApplyTest2 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscApplyTest8 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 8)
endif ()

## composites ##

waLBerla_compile_test(FILES composites/ConcatenatedOperator.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ConcatenatedOperator)

waLBerla_compile_test(FILES composites/P0P1HelperFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest1 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 1)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest3 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 3)

if (HYTEG_BUILD_WITH_PETSC)
    #waLBerla_compile_test(FILES composites/P1P0Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME P1P0Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P0Stokes2DPetscSolveTest> PROCESSES 1)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1StokesPetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1StokesPetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1StokesPetscSymmetryTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P1StokesPetscSymmetryTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P2Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P2Stokes2DPetscSolveTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P2Stokes3DPetscSolveTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME P2P2Stokes3DPetscSolveTest)
endif ()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES convergence/P2P1Stokes2DGKBConvergenceTest.cpp DEPENDS hyteg core )
    waLBerla_execute_test(NAME P2P1Stokes2DGKBConvergenceTest)
endif()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES convergence/P2P1ChannelGKBConvergenceTest.cpp DEPENDS hyteg core )
    waLBerla_execute_test(NAME P2P1ChannelGKBConvergenceTest)
endif()

waLBerla_compile_test(FILES composites/P2P1ElementwiseStokesOperatorTest.cpp DEPENDS hyteg core )
waLBerla_execute_test(NAME P2P1ElementwiseStokesOperatorTest)

############################################

## VectorFunctions ##
waLBerla_compile_test(FILES VectorFunctions/VectorFunctionBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VectorFunctionBasicTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES VectorFunctions/VectorFunctionPetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME VectorFunctionPetscTest)
endif ()

waLBerla_compile_test(FILES operators/VectorOperatorBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VectorOperatorBasicTest)

## Other ##


waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME vCycleTest1 COMMAND $<TARGET_FILE:vCycleTest>)



waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES LevelinfoTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME LevelinfoTest)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageTest1 COMMAND $<TARGET_FILE:PrimitiveStorageTest>)
waLBerla_execute_test(NAME PrimitiveStorageTest3 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageTest8 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 8)

waLBerla_compile_test(FILES PrimitiveRelativeOrderingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveRelativeOrderingTest)

waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BufferedCommunicationTest1 COMMAND $<TARGET_FILE:BufferedCommunicationTest>)
waLBerla_execute_test(NAME BufferedCommunicationTest3 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 3)
waLBerla_execute_test(NAME BufferedCommunicationTest8 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveMigrationTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationTest>)
waLBerla_execute_test(NAME PrimitiveMigrationTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationTest8 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/ParallelRoundRobinTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ParallelRoundRobinTest1 COMMAND $<TARGET_FILE:ParallelRoundRobinTest>)
waLBerla_execute_test(NAME ParallelRoundRobinTest3 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 3)
waLBerla_execute_test(NAME ParallelRoundRobinTest8 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationMatMulTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest>)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 4)

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest>)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 4)

waLBerla_compile_test(FILES adaptivity/PrimitiveStorageCopyTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest1 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest>)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest3 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest8 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 8)

waLBerla_compile_test(FILES adaptivity/AgglomerationMatMulTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME AgglomerationMatMulTest1 COMMAND $<TARGET_FILE:AgglomerationMatMulTest>)
waLBerla_execute_test(NAME AgglomerationMatMulTest2 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationMatMulTest3 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationMatMulTest4 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 4 )

waLBerla_compile_test(FILES adaptivity/AgglomerationConvergenceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME AgglomerationConvergenceTest1 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest>)
waLBerla_execute_test(NAME AgglomerationConvergenceTest2 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationConvergenceTest3 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationConvergenceTest4 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 4 )

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES adaptivity/AgglomerationSupermanFactorizationTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME AgglomerationSupermanFactorizationTest4 COMMAND $<TARGET_FILE:AgglomerationSupermanFactorizationTest> PROCESSES 4 )
endif ()

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

waLBerla_compile_test(FILES functions/FunctionMemorySerializationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest>)
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8)

waLBerla_compile_test(FILES functions/FunctionMultiStoreTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMultiStoreTest)

waLBerla_compile_test(FILES functions/FunctionRuleOfThreeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest1 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest>)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest2 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest8 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 8)

## Blending ##
waLBerla_compile_test(FILES blending/AffineMap3DTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME AffineMap3DTest1 COMMAND $<TARGET_FILE:AffineMap3DTest>)
endif ()


waLBerla_compile_test(FILES blending/AnnulusMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME AnnulusMapTest1 COMMAND $<TARGET_FILE:AnnulusMapTest>)

waLBerla_compile_test(FILES blending/IcosahedralShellMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IcosahedralShellMapTest COMMAND $<TARGET_FILE:IcosahedralShellMapTest>)

waLBerla_compile_test(FILES blending/InverseBlending+PrimitiveIDFindingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME InverseBlending+PrimitiveIDFindingTest)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (NOT WALBERLA_CXX_COMPILER_IS_INTEL)
        waLBerla_compile_test(FILES blending/P2P1StokesElementwiseAnnulusBlendingTest.cpp DEPENDS hyteg core)
        waLBerla_execute_test(NAME P2P1StokesElementwiseAnnulusBlendingTest COMMAND $<TARGET_FILE:P2P1StokesElementwiseAnnulusBlendingTest>)
    endif ()
endif ()

waLBerla_compile_test(FILES blending/PolygonDistanceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PolygonDistanceTest)

waLBerla_compile_test(FILES blending/TokamakMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME TokamakMapTest)

## Elementwise Operators ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES operators/ElementwiseOperatorPetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest>)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest2 COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest> PROCESSES 2)
endif ()

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

## numeric tools ##
waLBerla_compile_test(FILES numerictools/SpectrumEstimationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SpectrumEstimationTest)
waLBerla_execute_test(NAME SpectrumEstimationTestMPI COMMAND $<TARGET_FILE:SpectrumEstimationTest> PROCESSES 3)

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

## Free slip ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES freeslip/FreeslipPetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME FreeslipPetscApplyTest)
endif ()

waLBerla_compile_test(FILES freeslip/FreeslipRectangularChannelTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FreeslipRectangularChannelTest)
waLBerla_execute_test(NAME FreeslipRectangularChannelTestMPI COMMAND $<TARGET_FILE:FreeslipRectangularChannelTest> PROCESSES 2)

waLBerla_compile_test(FILES freeslip/ProjectNormalTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ProjectNormalTest)
waLBerla_execute_test(NAME ProjectNormalTestMPI COMMAND $<TARGET_FILE:ProjectNormalTest> PROCESSES 2)

## Trilinos ##
if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES trilinos/TrilinosVersionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosVersionTest)

    waLBerla_compile_test(FILES trilinos/TrilinosVectorTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosVectorTest)

    waLBerla_compile_test(FILES trilinos/TrilinosSparseMatrixTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosSparseMatrixTest)

    waLBerla_compile_test(FILES trilinos/TrilinosDirectSolverTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosDirectSolverTest)
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

waLBerla_compile_test(FILES operators/OperatorWrapperBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME OperatorWrapperBasicTest)

waLBerla_compile_test(FILES operators/BlockOperatorBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BlockOperatorBasicTest)

waLBerla_compile_test(FILES operators/BlockOperatorGSSmoothableTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BlockOperatorGSSmoothableTest)

## PETSc & Matrix Assembly ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Petsc/PetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest>)
    waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2)

    waLBerla_compile_test(FILES Petsc/PetscMatrixAssemblyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest>)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest2 COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest> PROCESSES 2)

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
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 2)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 4)

    waLBerla_compile_test(FILES dg/DGPoisson3DConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPoisson3DConvergenceTest)
    # Comment-in once macro-cell to macro-cell packing is implemented
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 2)
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 4)

    waLBerla_compile_test(FILES dg/DGLaplaceSymmetryTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGLaplaceSymmetryTest)

    waLBerla_compile_test(FILES dg/DGInterpolateEvaluateTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGInterpolateEvaluateTest)

    #waLBerla_compile_test(FILES egfunctionspace/EGFunctionspaceTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME EGFunctionspaceTest)

   # waLBerla_compile_test(FILES egfunctionspace/EGSymmetryTest.cpp DEPENDS hyteg core)
   # waLBerla_execute_test(NAME EGSymmetryTest)

   # waLBerla_compile_test(FILES egfunctionspace/EGVectorLaplaceConvergenceTests.cpp DEPENDS hyteg core)
   # waLBerla_execute_test(NAME EGVectorLaplaceConvergenceTests)

    #waLBerla_compile_test(FILES egfunctionspace/EGStokesConvergenceTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME EGStokesConvergenceTest COMMAND $<TARGET_FILE:EGStokesConvergenceTest> -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space)

    #waLBerla_compile_test(FILES egfunctionspace/EGNietscheBCTests.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME EGNietscheBCTests )

    #waLBerla_compile_test(FILES egfunctionspace/EGEpsilonConvergenceTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME EGEpsilonConvergenceTest COMMAND $<TARGET_FILE:EGEpsilonConvergenceTest> -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space)

    #waLBerla_compile_test(FILES egfunctionspace/EGFormTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME EGFormTest)

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

## Basic Smoothing and Solving
waLBerla_compile_test(FILES smooth+solve/SmoothableTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SmoothableTest)

waLBerla_compile_test(FILES smooth+solve/VectorToVectorOperatorCGTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VectorToVectorOperatorCGTest)

waLBerla_compile_test(FILES smooth+solve/VectorToVectorOperatorChebyshevTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VectorToVectorOperatorChebyshevTest)

## Structural
waLBerla_compile_test(FILES PrintInfoTest.cpp)
waLBerla_execute_test(NAME PrintInfoTest)

## adaptive refinement w/ hanging nodes
waLBerla_compile_test(FILES adaptiverefinementhangingnodes/ARHangingNodesStorageTest.cpp)
waLBerla_execute_test(NAME ARHangingNodesStorageTest)


if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES adaptiverefinementhangingnodes/P1AdaptiveRefinementEqualityTest.cpp)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME P1AdaptiveRefinementEqualityTest)
    endif()

    waLBerla_compile_test(FILES adaptiverefinementhangingnodes/DGAdaptiveRefinementEqualityTest.cpp)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME DGAdaptiveRefinementEqualityTest)
    endif()
endif ()

waLBerla_compile_test(FILES adaptiverefinementhangingnodes/MacroVolumeRefinementMicroIdxConversionTest.cpp)
waLBerla_execute_test(NAME MacroVolumeRefinementMicroIdxConversionTest)

waLBerla_compile_test(FILES eigen/EigenBasicTest.cpp)
waLBerla_execute_test(NAME EigenBasicTest)
waLBerla_execute_test(NAME EigenBasicTestMPI COMMAND $<TARGET_FILE:EigenBasicTest> PROCESSES 2)

## N1E1 ##
# For some reason executing some tests takes forever with clang and gcc in Debug mode.
set(N1E1_CONDITIONAL_LONGRUN "")
if (     CMAKE_BUILD_TYPE MATCHES "Debug"
     AND (WALBERLA_CXX_COMPILER_IS_CLANG OR WALBERLA_CXX_COMPILER_IS_GNU)
   )
    set(N1E1_CONDITIONAL_LONGRUN "longrun")
endif ()

waLBerla_compile_test(FILES N1E1/N1E1AddTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AddTest)

waLBerla_compile_test(FILES N1E1/N1E1AssignTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AssignTest)

waLBerla_compile_test(FILES N1E1/N1E1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1EvaluateTest)

waLBerla_compile_test(FILES N1E1/N1E1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1InterpolateTest)

waLBerla_compile_test(FILES N1E1/N1E1AdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AdditiveCommunicationTest)

waLBerla_compile_test(FILES N1E1/N1E1CurlCurlNullSpaceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlNullSpaceTest)

waLBerla_compile_test(FILES N1E1/N1E1MassTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1MassTest)

waLBerla_compile_test(FILES N1E1/N1E1FormCurlCurlTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1FormCurlCurlTest)

waLBerla_compile_test(FILES N1E1/N1E1CurlCurlConvergenceTest.cpp N1E1/common.hpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1CurlCurlConvergenceTest> LABELS longrun PROCESSES 4)

waLBerla_compile_test(FILES N1E1/N1E1ProlongationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1ProlongationTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1RestrictionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1RestrictionTest LABELS ${N1E1_CONDITIONAL_LONGRUN})
endif ()

waLBerla_compile_test(FILES N1E1/P1toN1E1GradientTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1toN1E1GradientTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1toP1LiftingTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1toP1LiftingTest LABELS ${N1E1_CONDITIONAL_LONGRUN})
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1/N1E1InverseDiagonalTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1InverseDiagonalTest)
endif ()

waLBerla_compile_test(FILES N1E1/HybridSmootherTest.cpp N1E1/common.hpp DEPENDS hyteg core)
waLBerla_execute_test(NAME HybridSmootherTest LABELS ${N1E1_CONDITIONAL_LONGRUN})

waLBerla_compile_test(FILES N1E1/N1E1GMGhIndependenceTest.cpp N1E1/common.hpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1GMGhIndependenceTest LABELS longrun)

waLBerla_compile_test(FILES N1E1/N1E1GMGConvergenceTest.cpp N1E1/common.hpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1GMGConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1GMGConvergenceTest> LABELS longrun PROCESSES 4)

## types ##
waLBerla_compile_test(FILES types/MatrixConstructFromEigenTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MatrixConstructFromEigenTest)

waLBerla_compile_test(FILES mixedPrecision/basicMixedPrecisionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME basicMixedPrecisionTest)

waLBerla_compile_test(FILES mixedPrecision/mixedPrecisionIterativeRefinement.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME mixedPrecisionIterativeRefinement)
