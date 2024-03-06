waLBerla_compile_test(FILES convergence/P1VectorLaplaceCGTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P1VectorLaplaceCGTest)

waLBerla_compile_test(FILES convergence/P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P1StokesMinResConvergenceTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES convergence/P1P1StokesUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1StokesUzawaConvergenceTest)
endif ()

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_compile_test(FILES convergence/P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1Stokes3DMinResConvergenceTest)

    waLBerla_compile_test(FILES convergence/P2P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes3DMinResConvergenceTest)

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES convergence/P1P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
        waLBerla_execute_test(NAME P1P1Stokes3DUzawaConvergenceTest)
    endif ()

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES convergence/P2P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
        waLBerla_execute_test(NAME P2P1Stokes3DUzawaConvergenceTest)
    endif ()
endif ()

waLBerla_compile_test(FILES convergence/P2ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P2ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES convergence/ElementwiseEpsilonMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-stokes)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTest)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTestLongrun COMMAND $<TARGET_FILE:ElementwiseEpsilonMinResConvergenceTest> --longrun LABELS longrun)
    endif ()
endif ()

waLBerla_compile_test(FILES convergence/P2P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1StokesMinResConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1StokesSchurCGConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1StokesSchurCGConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1ElementwiseUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1ElementwiseUzawaConvergenceTest)

waLBerla_compile_test(FILES convergence/P2P1UzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1UzawaConvergenceTest)

waLBerla_compile_test(FILES convergence/StokesVarViscPopovTest.cpp DEPENDS hyteg core mixed_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME StokesVarViscPopovTest)
    endif ()
endif ()

waLBerla_compile_test(FILES dataexport+import/VTKOutputTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VTKOutputTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Indexing/CouplingCountTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME CouplingCountTest COMMAND $<TARGET_FILE:CouplingCountTest>)
    waLBerla_execute_test(NAME ParallelCouplingCountTest3 COMMAND $<TARGET_FILE:CouplingCountTest> PROCESSES 3)
endif ()

waLBerla_compile_test(FILES operators/DivergenceOperatorTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-divergence)
waLBerla_execute_test(NAME DivergenceOperatorTest)

waLBerla_compile_test(FILES operators/ElementwiseOperatorAdditiveApplyTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME ElementwiseOperatorAdditiveApplyTest)

waLBerla_compile_test(FILES operators/ViscousOperatorsTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-viscousblock)
waLBerla_execute_test(NAME ViscousOperatorsTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1StokesPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1StokesPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P1P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P1StokesPetscSymmetryTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1StokesPetscSymmetryTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P2Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P2Stokes2DPetscSolveTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES composites/P2P2Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P2Stokes3DPetscSolveTest)
endif ()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES convergence/P2P1Stokes2DGKBConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes2DGKBConvergenceTest)
endif()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES convergence/P2P1ChannelGKBConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1ChannelGKBConvergenceTest)
endif()

waLBerla_compile_test(FILES composites/P2P1ElementwiseStokesOperatorTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-stokes)
waLBerla_execute_test(NAME P2P1ElementwiseStokesOperatorTest)

waLBerla_compile_test(FILES operators/VectorOperatorBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorOperatorBasicTest)

#if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (NOT WALBERLA_CXX_COMPILER_IS_INTEL)
        waLBerla_compile_test(FILES blending/P2P1StokesElementwiseAnnulusBlendingTest.cpp DEPENDS hyteg core mixed_operator)
        waLBerla_execute_test(NAME P2P1StokesElementwiseAnnulusBlendingTest COMMAND $<TARGET_FILE:P2P1StokesElementwiseAnnulusBlendingTest>)
    endif ()
#endif ()

## Elementwise Operators ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES operators/ElementwiseOperatorPetscTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest>)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest2 COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest> PROCESSES 2)
endif ()

## Free slip ##
if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES freeslip/FreeslipPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME FreeslipPetscApplyTest)
endif ()

waLBerla_compile_test(FILES freeslip/FreeslipRectangularChannelTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME FreeslipRectangularChannelTest)
waLBerla_execute_test(NAME FreeslipRectangularChannelTestMPI COMMAND $<TARGET_FILE:FreeslipRectangularChannelTest> PROCESSES 2)

if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES trilinos/TrilinosSparseMatrixTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME TrilinosSparseMatrixTest)

    waLBerla_compile_test(FILES trilinos/TrilinosDirectSolverTest.cpp DEPENDS hyteg core mixed_operator)
    #TODO: renable test one the spack installed trilinos contains MUMPS (spack v0.20)
    #waLBerla_execute_test(NAME TrilinosDirectSolverTest)
endif ()

waLBerla_compile_test(FILES operators/OperatorWrapperBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME OperatorWrapperBasicTest)

waLBerla_compile_test(FILES operators/BlockOperatorBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME BlockOperatorBasicTest)

waLBerla_compile_test(FILES operators/BlockOperatorGSSmoothableTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME BlockOperatorGSSmoothableTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES Petsc/PetscMatrixAssemblyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest>)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest2 COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest> PROCESSES 2)

    waLBerla_compile_test(FILES egfunctionspace/EGSymmetryTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGSymmetryTest)

    waLBerla_compile_test(FILES egfunctionspace/EGApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGApplyTest)

    waLBerla_compile_test(FILES egfunctionspace/EGMatfreeConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGMatfreeConvergenceTest COMMAND $<TARGET_FILE:EGMatfreeConvergenceTest> -ksp_type
            minres
            -pc_type
            none
            -ksp_monitor)

    waLBerla_compile_test(FILES egfunctionspace/EGStokesConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGStokesConvergenceTest COMMAND $<TARGET_FILE:EGStokesConvergenceTest> )

    waLBerla_compile_test(FILES egfunctionspace/EGEpsilonConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGEpsilonConvergenceTest COMMAND $<TARGET_FILE:EGEpsilonConvergenceTest> )
endif ()

## Basic Smoothing and Solving
waLBerla_compile_test(FILES smooth+solve/SmoothableTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME SmoothableTest)

waLBerla_compile_test(FILES smooth+solve/VectorToVectorOperatorCGTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorToVectorOperatorCGTest)

waLBerla_compile_test(FILES smooth+solve/VectorToVectorOperatorChebyshevTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorToVectorOperatorChebyshevTest)

waLBerla_compile_test(FILES N1E1/HybridSmootherTest.cpp ../hyteg/N1E1/common.hpp ../hyteg/N1E1/common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME HybridSmootherTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})

waLBerla_compile_test(FILES N1E1/N1E1GMGhIndependenceTest.cpp ../hyteg/N1E1/common.hpp ../hyteg/N1E1/common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME N1E1GMGhIndependenceTest LABELS longrun LABELS n1e1)

waLBerla_compile_test(FILES N1E1/N1E1GMGConvergenceTest.cpp ../hyteg/N1E1/common.hpp ../hyteg/N1E1/common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME N1E1GMGConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1GMGConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

