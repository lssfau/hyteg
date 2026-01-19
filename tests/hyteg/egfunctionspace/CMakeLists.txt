if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( EGFunctionspaceTest EGFunctionspaceTest.cpp )
    target_link_libraries       ( EGFunctionspaceTest hyteg walberla::core )
    waLBerla_execute_test(NAME EGFunctionspaceTest)

    waLBerla_add_test_executable( EGSymmetryTest EGSymmetryTest.cpp )
    target_link_libraries       ( EGSymmetryTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME EGSymmetryTest)

    waLBerla_add_test_executable( EGApplyTest EGApplyTest.cpp )
    target_link_libraries       ( EGApplyTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME EGApplyTest)

    waLBerla_add_test_executable( EGMatfreeConvergenceTest EGMatfreeConvergenceTest.cpp )
    target_link_libraries       ( EGMatfreeConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME EGMatfreeConvergenceTest COMMAND $<TARGET_FILE:EGMatfreeConvergenceTest> -ksp_type
          minres
          -pc_type
          none
          -ksp_monitor)

    waLBerla_add_test_executable( EGStokesConvergenceTest EGStokesConvergenceTest.cpp )
    target_link_libraries       ( EGStokesConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME EGStokesConvergenceTest COMMAND $<TARGET_FILE:EGStokesConvergenceTest> )

    waLBerla_add_test_executable( EGEpsilonConvergenceTest EGEpsilonConvergenceTest.cpp )
    target_link_libraries       ( EGEpsilonConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME EGEpsilonConvergenceTest COMMAND $<TARGET_FILE:EGEpsilonConvergenceTest> )
endif ()