waLBerla_compile_test(FILES P0P1HelperFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest1 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 1)
waLBerla_execute_test(NAME P0P1HelperFunctionsTest3 COMMAND $<TARGET_FILE:P0P1HelperFunctionsTest> PROCESSES 3)

if (HYTEG_BUILD_WITH_PETSC)
    #waLBerla_compile_test(FILES P1P0Stokes2DPetscSolveTest.cpp DEPENDS hyteg core)
    #waLBerla_execute_test(NAME P1P0Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P0Stokes2DPetscSolveTest> PROCESSES 1)
endif ()

waLBerla_compile_test(FILES ConcatenatedOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME ConcatenatedOperatorTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1P1StokesPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P1P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P1StokesPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest1 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1StokesPetscApplyTest2 COMMAND $<TARGET_FILE:P2P1StokesPetscApplyTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P1Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes2DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes2DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P1P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P1Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest1 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2P1Stokes3DPetscSolveTest2 COMMAND $<TARGET_FILE:P2P1Stokes3DPetscSolveTest> PROCESSES 2)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P1StokesPetscSymmetryTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1StokesPetscSymmetryTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P2Stokes2DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P2Stokes2DPetscSolveTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2P2Stokes3DPetscSolveTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P2Stokes3DPetscSolveTest)
endif ()

waLBerla_compile_test(FILES P2P1ElementwiseStokesOperatorTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-stokes)
waLBerla_execute_test(NAME P2P1ElementwiseStokesOperatorTest)

