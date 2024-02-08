## Convergence ##
waLBerla_compile_test(FILES convergence/P1CGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1CGConvergenceTest)
waLBerla_execute_test(NAME P1CGConvergenceTestMPI COMMAND $<TARGET_FILE:P1CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1CG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1CG3DConvergenceTest)
waLBerla_execute_test(NAME P1CG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1CG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1GMGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMGConvergenceTest)
waLBerla_execute_test(NAME P1GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1FASConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1FASConvergenceTest)
waLBerla_execute_test(NAME P1FASConvergenceTestMPI COMMAND $<TARGET_FILE:P1FASConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1TauExtrapolationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTest)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTestMPI COMMAND $<TARGET_FILE:P1TauExtrapolationConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1GMG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMG3DConvergenceTest)
waLBerla_execute_test(NAME P1GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMGConvergenceTest)
waLBerla_execute_test(NAME P2GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMG3DConvergenceTest)
waLBerla_execute_test(NAME P2GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMG3DConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P1GSConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1GSConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P1JacobiConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1JacobiConvergenceTest)

waLBerla_compile_test(FILES convergence/P1MinResConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1MinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P2CGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2CGConvergenceTest)
waLBerla_execute_test(NAME P2CGConvergenceTestMPI COMMAND $<TARGET_FILE:P2CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P1ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P1ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2P1FGMRESStokesTest.cpp DEPENDS hyteg core constant_stencil_operator mixed_operator)
waLBerla_execute_test(NAME P2P1FGMRESStokesTest)

waLBerla_compile_test(FILES convergence/P2P1GMRESStokesTest.cpp DEPENDS hyteg core constant_stencil_operator mixed_operator)
waLBerla_execute_test(NAME P2P1GMRESStokesTest)

waLBerla_compile_test(FILES convergence/P2GMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMRESConvergenceTest)

waLBerla_compile_test(FILES convergence/P2FGMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2FGMRESConvergenceTest)

waLBerla_compile_test(FILES convergence/P1GMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMRESConvergenceTest)

waLBerla_compile_test(FILES convergence/P1FGMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1FGMRESConvergenceTest)

waLBerla_compile_test(FILES convergence/P2SmoothConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2SmoothConvergenceTest)

waLBerla_compile_test(FILES convergence/P23DSmoothConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P23DSmoothConvergenceTest)

waLBerla_compile_test(FILES convergence/P2JacobiConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2JacobiConvergenceTest)

waLBerla_compile_test(FILES convergence/P2MinResConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2MinResConvergenceTest)
waLBerla_execute_test(NAME P2MinResConvergenceTestMPI COMMAND $<TARGET_FILE:P2MinResConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGHRefinementConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2GMGHRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGHRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGHRefinementConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/P2GMGPRefinementConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2GMGPRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGPRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGPRefinementConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P2GSConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2GSConvergenceTest)
endif()

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES convergence/P1ChebyshevSmootherConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1ChebyshevSmootherConvergenceTest)
endif()

waLBerla_compile_test(FILES convergence/P2UnsteadyDiffusion2DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2UnsteadyDiffusion2DTest)

waLBerla_compile_test(FILES convergence/P2UnsteadyDiffusionBlending2DTest.cpp DEPENDS hyteg core constant_stencil_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME P2UnsteadyDiffusionBlending2DTest)
    endif ()
endif ()

waLBerla_compile_test(FILES convection/Upwind2DCircularConvectionTest.cpp DEPENDS hyteg core)
waLBerla_compile_test(FILES convection/AlgebraicUpwind3DCircularConvectionTest.cpp DEPENDS hyteg core constant_stencil_operator)

waLBerla_compile_test(FILES forms/P2LinearCombinationFormTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2LinearCombinationFormTest)

waLBerla_compile_test(FILES forms/RowSumTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME RowSumTest)

## meshes ##
waLBerla_compile_test(FILES meshes/RefineCoarseMeshCGTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME RefineCoarseMeshCGTest)

waLBerla_compile_test(FILES vertexdofspace/VertexDoFFunctionTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME VertexDoFFunctionTest)

waLBerla_compile_test(FILES operators/EdgeDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME EdgeDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES operators/MassOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME MassOperatorTest)
    waLBerla_execute_test(NAME MassOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 2)
endif ()

waLBerla_compile_test(FILES operators/DiagonalNonConstantOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME DiagonalNonConstantOperatorTest)

waLBerla_compile_test(FILES P1/P1LaplaceOperator3DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1/P1PetscApplyTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1PetscApplyTest1 COMMAND $<TARGET_FILE:P1PetscApplyTest>)
    waLBerla_execute_test(NAME P1PetscApplyTest2 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscApplyTest8 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1/P1PetscSolveTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1PetscSolveTest1 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1PetscSolveTest2 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscSolveTest3 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 8)
endif ()

waLBerla_compile_test(FILES P2/P2GSTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GSTest)

waLBerla_compile_test(FILES P2/P2SimpleSmoothTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2SimpleSmoothTest)

waLBerla_compile_test(FILES P2/P2JacobiGSSmoothTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2JacobiGSSmoothTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2LaplacePetscSymmetryTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2LaplacePetscSymmetryTest)
endif ()

waLBerla_compile_test(FILES P2/P2LaplaceOperator3DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscSolveTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscSolveTest1 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2PetscSolveTest2 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 2)
    if (NOT WALBERLA_CXX_COMPILER_IS_INTELLLVM)
        waLBerla_execute_test(NAME P2PetscSolveTest3 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 8)
    endif ()
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscTest1 COMMAND $<TARGET_FILE:P2PetscTest>)
    waLBerla_execute_test(NAME P2PetscTest2 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscTest8 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2/P2PetscApplyTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscApplyTest1 COMMAND $<TARGET_FILE:P2PetscApplyTest>)
    waLBerla_execute_test(NAME P2PetscApplyTest2 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscApplyTest8 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 8)
endif ()

waLBerla_compile_test(FILES composites/ConcatenatedOperator.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME ConcatenatedOperator)

## VectorFunctions ##
waLBerla_compile_test(FILES VectorFunctions/VectorFunctionBasicTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME VectorFunctionBasicTest)

waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME vCycleTest1 COMMAND $<TARGET_FILE:vCycleTest>)

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationMatMulTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest>)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 4)

waLBerla_compile_test(FILES adaptivity/PrimitiveMigrationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest>)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 4)

waLBerla_compile_test(FILES adaptivity/AgglomerationMatMulTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME AgglomerationMatMulTest1 COMMAND $<TARGET_FILE:AgglomerationMatMulTest>)
waLBerla_execute_test(NAME AgglomerationMatMulTest2 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationMatMulTest3 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationMatMulTest4 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 4 )

waLBerla_compile_test(FILES adaptivity/AgglomerationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME AgglomerationConvergenceTest1 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest>)
waLBerla_execute_test(NAME AgglomerationConvergenceTest2 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationConvergenceTest3 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationConvergenceTest4 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 4 )

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES adaptivity/AgglomerationSupermanFactorizationTest.cpp DEPENDS hyteg core constant_stencil_operator)
    #waLBerla_execute_test(NAME AgglomerationSupermanFactorizationTest4 COMMAND $<TARGET_FILE:AgglomerationSupermanFactorizationTest> PROCESSES 4 )
endif ()

waLBerla_compile_test(FILES functions/FunctionMemorySerializationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest>)
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8)

## numeric tools ##
waLBerla_compile_test(FILES numerictools/SpectrumEstimationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME SpectrumEstimationTest)
waLBerla_execute_test(NAME SpectrumEstimationTestMPI COMMAND $<TARGET_FILE:SpectrumEstimationTest> PROCESSES 3)

if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES trilinos/TrilinosVectorTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME TrilinosVectorTest)
endif ()

## PETSc & Matrix Assembly ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Petsc/PetscTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest>)
    waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2)

    waLBerla_compile_test(FILES Petsc/PetscApplyDirichletBoundaryUIDTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME PetscApplyDirichletBoundaryUIDTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES adaptiverefinementhangingnodes/P1AdaptiveRefinementEqualityTest.cpp DEPENDS hyteg constant_stencil_operator)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME P1AdaptiveRefinementEqualityTest)
    endif()

    waLBerla_compile_test(FILES adaptiverefinementhangingnodes/DGAdaptiveRefinementEqualityTest.cpp DEPENDS hyteg constant_stencil_operator)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME DGAdaptiveRefinementEqualityTest)
    endif()
endif ()

waLBerla_compile_test(FILES N1E1/N1E1BasisTransformationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME N1E1BasisTransformationTest LABELS n1e1)

waLBerla_compile_test(FILES mixedPrecision/mixedPrecisionIterativeRefinement.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME mixedPrecisionIterativeRefinement)

