#if( HYTEG_BUILD_WITH_PETSC )
#    waLBerla_compile_test(FILES P1WithDGFromPoisson2DConvergenceTest.cpp DEPENDS hyteg core)
#    waLBerla_execute_test(NAME P1WithDGFromPoisson2DConvergenceTest)
#endif()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES DGPoisson2DConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 2 LABELS longrun)
    waLBerla_execute_test(NAME DGPoisson2DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson2DConvergenceTest> PROCESSES 4 LABELS longrun)

    waLBerla_compile_test(FILES DGPoisson3DConvergenceTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPoisson3DConvergenceTest)
    # Comment-in once macro-cell to macro-cell packing is implemented
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest2 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 2)
    # waLBerla_execute_test(NAME DGPoisson3DConvergenceTest4 COMMAND $<TARGET_FILE:DGPoisson3DConvergenceTest> PROCESSES 4)

    waLBerla_compile_test(FILES DGLaplaceSymmetryTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGLaplaceSymmetryTest)

    waLBerla_compile_test(FILES DGPetscApplyTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGPetscApplyTest)

    waLBerla_compile_test(FILES DGGlobalRefinementGridConvTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGGlobalRefinementGridConvTest)

    waLBerla_compile_test(FILES DGAdaptiveRefinementInterpolateEvaluateTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGAdaptiveRefinementInterpolateEvaluateTest)

    waLBerla_compile_test(FILES DGAdaptiveRefinementGridConvTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGAdaptiveRefinementGridConvTest)

    waLBerla_compile_test(FILES DGProlongationTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME DGProlongationTest)
endif ()

waLBerla_compile_test(FILES DGInterpolateEvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME DGInterpolateEvaluateTest)
