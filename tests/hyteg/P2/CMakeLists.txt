waLBerla_add_test_executable( P2P1TransferTest P2P1TransferTest.cpp )
target_link_libraries       ( P2P1TransferTest hyteg walberla::core )
waLBerla_execute_test(NAME P2P1TransferTest)

waLBerla_add_test_executable( P2VertexInterpolateTest P2VertexInterpolateTest.cpp )
target_link_libraries       ( P2VertexInterpolateTest hyteg walberla::core )
waLBerla_execute_test(NAME P2VertexInterpolateTest)

waLBerla_add_test_executable( P2RestrictTest P2RestrictTest.cpp )
target_link_libraries       ( P2RestrictTest hyteg walberla::core )
waLBerla_execute_test(NAME P2RestrictTest)

waLBerla_add_test_executable( P2ProlongateTest P2ProlongateTest.cpp )
target_link_libraries       ( P2ProlongateTest hyteg walberla::core )
waLBerla_execute_test(NAME P2ProlongateTest)

waLBerla_add_test_executable( P2QuadraticRestriction3DTest P2QuadraticRestriction3DTest.cpp )
target_link_libraries       ( P2QuadraticRestriction3DTest hyteg walberla::core )
waLBerla_execute_test(NAME P2QuadraticRestriction3DTest)

waLBerla_add_test_executable( P2QuadraticProlongation3DTest P2QuadraticProlongation3DTest.cpp )
target_link_libraries       ( P2QuadraticProlongation3DTest hyteg walberla::core )
waLBerla_execute_test(NAME P2QuadraticProlongation3DTest)

waLBerla_add_test_executable( P2BasicFunctionsTest P2BasicFunctionsTest.cpp )
target_link_libraries       ( P2BasicFunctionsTest hyteg walberla::core )
waLBerla_execute_test(NAME P2BasicFunctionsTest)

waLBerla_add_test_executable( P2SwapTest P2SwapTest.cpp )
target_link_libraries       ( P2SwapTest hyteg walberla::core )
waLBerla_execute_test(NAME P2SwapTest)

waLBerla_add_test_executable( P2FunctionTest P2FunctionTest.cpp )
target_link_libraries       ( P2FunctionTest hyteg walberla::core )
waLBerla_execute_test(NAME P2FunctionTest)

target_sources( hyteg
    PRIVATE
    P2PlusBubbleElementwiseMass_AffineMap2D_float64.cpp
    P2PlusBubbleElementwiseMass_AffineMap2D_float64.hpp
)

waLBerla_add_test_executable( P2PlusBubbleFunctionTest P2PlusBubbleFunctionTest.cpp )
target_link_libraries       ( P2PlusBubbleFunctionTest hyteg opgen-diffusion opgen-mass walberla::core )
waLBerla_execute_test(NAME P2PlusBubbleFunctionTest1 COMMAND $<TARGET_FILE:P2PlusBubbleFunctionTest> PROCESSES 1)
waLBerla_execute_test(NAME P2PlusBubbleFunctionTest2 COMMAND $<TARGET_FILE:P2PlusBubbleFunctionTest> PROCESSES 2)

waLBerla_add_test_executable( P2Level0Test P2Level0Test.cpp )
target_link_libraries       ( P2Level0Test hyteg walberla::core )
waLBerla_execute_test(NAME P2Level0Test1 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 1)
waLBerla_execute_test(NAME P2Level0Test2 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 2)

waLBerla_add_test_executable( P2EvaluateTest P2EvaluateTest.cpp )
target_link_libraries       ( P2EvaluateTest hyteg walberla::core )
waLBerla_execute_test(NAME P2EvaluateTest)

waLBerla_add_test_executable( P2GSTest P2GSTest.cpp )
target_link_libraries       ( P2GSTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2GSTest)

waLBerla_add_test_executable( P2SimpleSmoothTest P2SimpleSmoothTest.cpp )
target_link_libraries       ( P2SimpleSmoothTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2SimpleSmoothTest)

waLBerla_add_test_executable( P2JacobiGSSmoothTest P2JacobiGSSmoothTest.cpp )
target_link_libraries       ( P2JacobiGSSmoothTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2JacobiGSSmoothTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P2LaplacePetscSymmetryTest P2LaplacePetscSymmetryTest.cpp )
    target_link_libraries       ( P2LaplacePetscSymmetryTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P2LaplacePetscSymmetryTest)
endif ()

waLBerla_add_test_executable( P2LaplaceOperator3DTest P2LaplaceOperator3DTest.cpp )
target_link_libraries       ( P2LaplaceOperator3DTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P2PetscSolveTest P2PetscSolveTest.cpp )
    target_link_libraries       ( P2PetscSolveTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P2PetscSolveTest1 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2PetscSolveTest2 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 2)
    if (NOT WALBERLA_CXX_COMPILER_IS_INTELLLVM)
        waLBerla_execute_test(NAME P2PetscSolveTest3 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 8)
    endif ()
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P2PetscTest P2PetscTest.cpp )
    target_link_libraries       ( P2PetscTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P2PetscTest1 COMMAND $<TARGET_FILE:P2PetscTest>)
    waLBerla_execute_test(NAME P2PetscTest2 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscTest8 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P2PetscApplyTest P2PetscApplyTest.cpp )
    target_link_libraries       ( P2PetscApplyTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P2PetscApplyTest1 COMMAND $<TARGET_FILE:P2PetscApplyTest>)
    waLBerla_execute_test(NAME P2PetscApplyTest2 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscApplyTest8 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 8)
endif ()

