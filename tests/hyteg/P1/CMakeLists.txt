waLBerla_add_test_executable( P1InterpolateTest P1InterpolateTest.cpp )
target_link_libraries       ( P1InterpolateTest hyteg walberla::core )
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_add_test_executable( P1CommTest P1CommTest.cpp )
target_link_libraries       ( P1CommTest hyteg walberla::core )
waLBerla_execute_test(NAME P1CommTest)

waLBerla_add_test_executable( P1DataHandlingTest P1DataHandlingTest.cpp )
target_link_libraries       ( P1DataHandlingTest hyteg walberla::core )
waLBerla_execute_test(NAME P1DataHandlingTest)

waLBerla_add_test_executable( P1EvaluateTest P1EvaluateTest.cpp )
target_link_libraries       ( P1EvaluateTest hyteg walberla::core )
waLBerla_execute_test(NAME P1EvaluateTest)

waLBerla_add_test_executable( P1LaplaceOperator3DTest P1LaplaceOperator3DTest.cpp )
target_link_libraries       ( P1LaplaceOperator3DTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P1PetscApplyTest P1PetscApplyTest.cpp )
    target_link_libraries       ( P1PetscApplyTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P1PetscApplyTest1 COMMAND $<TARGET_FILE:P1PetscApplyTest>)
    waLBerla_execute_test(NAME P1PetscApplyTest2 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscApplyTest8 COMMAND $<TARGET_FILE:P1PetscApplyTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P1PetscSolveTest P1PetscSolveTest.cpp )
    target_link_libraries       ( P1PetscSolveTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P1PetscSolveTest1 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P1PetscSolveTest2 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 2)
    waLBerla_execute_test(NAME P1PetscSolveTest3 COMMAND $<TARGET_FILE:P1PetscSolveTest> PROCESSES 8)
endif ()

