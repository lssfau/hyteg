waLBerla_compile_test(FILES P2P1TransferTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2P1TransferTest)

waLBerla_compile_test(FILES P2VertexInterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2VertexInterpolateTest)

waLBerla_compile_test(FILES P2RestrictTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2RestrictTest)

waLBerla_compile_test(FILES P2ProlongateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2ProlongateTest)

waLBerla_compile_test(FILES P2QuadraticRestriction3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2QuadraticRestriction3DTest)

waLBerla_compile_test(FILES P2QuadraticProlongation3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2QuadraticProlongation3DTest)

waLBerla_compile_test(FILES P2BasicFunctionsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2BasicFunctionsTest)

waLBerla_compile_test(FILES P2SwapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2SwapTest)

waLBerla_compile_test(FILES P2FunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2FunctionTest)

target_sources( hyteg
    PRIVATE
    P2PlusBubbleElementwiseMass_AffineMap2D_float64.cpp
    P2PlusBubbleElementwiseMass_AffineMap2D_float64.hpp
)

waLBerla_compile_test(FILES P2PlusBubbleFunctionTest.cpp DEPENDS hyteg opgen-diffusion opgen-mass core)
waLBerla_execute_test(NAME P2PlusBubbleFunctionTest)

waLBerla_compile_test(FILES P2Level0Test.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2Level0Test1 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 1)
waLBerla_execute_test(NAME P2Level0Test2 COMMAND $<TARGET_FILE:P2Level0Test> PROCESSES 2)

waLBerla_compile_test(FILES P2EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P2EvaluateTest)

waLBerla_compile_test(FILES P2GSTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GSTest)

waLBerla_compile_test(FILES P2SimpleSmoothTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2SimpleSmoothTest)

waLBerla_compile_test(FILES P2JacobiGSSmoothTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2JacobiGSSmoothTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2LaplacePetscSymmetryTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2LaplacePetscSymmetryTest)
endif ()

waLBerla_compile_test(FILES P2LaplaceOperator3DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2LaplaceOperator3DTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2PetscSolveTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscSolveTest1 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 1)
    waLBerla_execute_test(NAME P2PetscSolveTest2 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 2)
    if (NOT WALBERLA_CXX_COMPILER_IS_INTELLLVM)
        waLBerla_execute_test(NAME P2PetscSolveTest3 COMMAND $<TARGET_FILE:P2PetscSolveTest> PROCESSES 8)
    endif ()
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2PetscTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscTest1 COMMAND $<TARGET_FILE:P2PetscTest>)
    waLBerla_execute_test(NAME P2PetscTest2 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscTest8 COMMAND $<TARGET_FILE:P2PetscTest> PROCESSES 8)
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P2PetscApplyTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2PetscApplyTest1 COMMAND $<TARGET_FILE:P2PetscApplyTest>)
    waLBerla_execute_test(NAME P2PetscApplyTest2 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 2)
    waLBerla_execute_test(NAME P2PetscApplyTest8 COMMAND $<TARGET_FILE:P2PetscApplyTest> PROCESSES 8)
endif ()

