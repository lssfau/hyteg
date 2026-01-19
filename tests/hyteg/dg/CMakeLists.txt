#if( HYTEG_BUILD_WITH_PETSC )
#    waLBerla_add_test_executable( P1WithDGFromPoisson2DConvergenceTest P1WithDGFromPoisson2DConvergenceTest.cpp )
#    target_link_libraries       ( P1WithDGFromPoisson2DConvergenceTest hyteg walberla::core )
#    waLBerla_execute_test(NAME P1WithDGFromPoisson2DConvergenceTest)
#endif()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( DGPoisson2DConvergenceTest DGPoisson2DConvergenceTest.cpp )
    target_link_libraries       ( DGPoisson2DConvergenceTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 2 LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 4 LABELS longrun)

    waLBerla_add_test_executable( DGPoisson3DConvergenceTest DGPoisson3DConvergenceTest.cpp )
    target_link_libraries       ( DGPoisson3DConvergenceTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGPoisson3DConvergenceTest)
    # Comment-in once macro-cell to macro-cell packing is implemented
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 2)
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 4)

    waLBerla_add_test_executable( DGLaplaceSymmetryTest DGLaplaceSymmetryTest.cpp )
    target_link_libraries       ( DGLaplaceSymmetryTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGLaplaceSymmetryTest)

    waLBerla_add_test_executable( DGPetscApplyTest DGPetscApplyTest.cpp )
    target_link_libraries       ( DGPetscApplyTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGPetscApplyTest)

    waLBerla_add_test_executable( DGGlobalRefinementGridConvTest DGGlobalRefinementGridConvTest.cpp )
    target_link_libraries       ( DGGlobalRefinementGridConvTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGGlobalRefinementGridConvTest)

    waLBerla_add_test_executable( DGAdaptiveRefinementInterpolateEvaluateTest DGAdaptiveRefinementInterpolateEvaluateTest.cpp )
    target_link_libraries       ( DGAdaptiveRefinementInterpolateEvaluateTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGAdaptiveRefinementInterpolateEvaluateTest)

    waLBerla_add_test_executable( DGAdaptiveRefinementGridConvTest DGAdaptiveRefinementGridConvTest.cpp )
    target_link_libraries       ( DGAdaptiveRefinementGridConvTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGAdaptiveRefinementGridConvTest)

    waLBerla_add_test_executable( DGProlongationTest DGProlongationTest.cpp )
    target_link_libraries       ( DGProlongationTest hyteg walberla::core )
    waLBerla_execute_test(NAME DGProlongationTest)
endif ()

waLBerla_add_test_executable( DGInterpolateEvaluateTest DGInterpolateEvaluateTest.cpp )
target_link_libraries       ( DGInterpolateEvaluateTest hyteg walberla::core )
waLBerla_execute_test(NAME DGInterpolateEvaluateTest)
