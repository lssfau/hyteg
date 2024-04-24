waLBerla_compile_test(FILES P1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES P1CommTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES P1DataHandlingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1DataHandlingTest)

waLBerla_compile_test(FILES P1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1EvaluateTest)

waLBerla_compile_test(FILES P1LaplaceOperator3DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1PetscApplyTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1PetscApplyTest1 COMMAND $<TARGET_FILE:P1PetscApplyTest>)
    waLBerla_execute_test(NAME P1PetscApplyTest2 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscApplyTest8 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1PetscSolveTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1PetscSolveTest1 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1PetscSolveTest2 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscSolveTest3 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 8)
endif ()

