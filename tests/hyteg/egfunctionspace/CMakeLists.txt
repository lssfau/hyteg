if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES EGFunctionspaceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME EGFunctionspaceTest)

    waLBerla_compile_test(FILES EGSymmetryTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGSymmetryTest)

    waLBerla_compile_test(FILES EGApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGApplyTest)

    waLBerla_compile_test(FILES EGMatfreeConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGMatfreeConvergenceTest COMMAND $<TARGET_FILE:EGMatfreeConvergenceTest> -ksp_type
          minres
          -pc_type
          none
          -ksp_monitor)

    waLBerla_compile_test(FILES EGStokesConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGStokesConvergenceTest COMMAND $<TARGET_FILE:EGStokesConvergenceTest> )

    waLBerla_compile_test(FILES EGEpsilonConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME EGEpsilonConvergenceTest COMMAND $<TARGET_FILE:EGEpsilonConvergenceTest> )
endif ()