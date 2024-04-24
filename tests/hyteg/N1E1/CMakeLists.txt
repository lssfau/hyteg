# For some reason executing some tests takes forever with clang and gcc in Debug mode.
set(N1E1_CONDITIONAL_LONGRUN "")
if (     CMAKE_BUILD_TYPE MATCHES "Debug"
      AND (WALBERLA_CXX_COMPILER_IS_CLANG OR WALBERLA_CXX_COMPILER_IS_GNU)
)
    set(N1E1_CONDITIONAL_LONGRUN "longrun")
endif ()

waLBerla_compile_test(FILES N1E1AddTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AddTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1AssignTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AssignTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1EvaluateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1EvaluateTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1InterpolateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1InterpolateTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1AdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1AdditiveCommunicationTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1CurlCurlNullSpaceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlNullSpaceTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1MassTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1MassTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1FormCurlCurlTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1FormCurlCurlTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1CurlCurlConvergenceTest.cpp common.hpp common.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1CurlCurlConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1CurlCurlConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

waLBerla_compile_test(FILES N1E1ProlongationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1ProlongationTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1RestrictionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1RestrictionTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

waLBerla_compile_test(FILES P1toN1E1GradientTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1toN1E1GradientTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1toP1LiftingTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1toP1LiftingTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES N1E1InverseDiagonalTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME N1E1InverseDiagonalTest LABELS n1e1)
endif ()

waLBerla_compile_test(FILES N1E1BlendingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1BlendingTest LABELS n1e1)

waLBerla_compile_test(FILES N1E1LinearFormBlendingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME N1E1LinearFormBlendingTest LABELS n1e1)

waLBerla_compile_test(FILES HybridSmootherTest.cpp common.hpp common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME HybridSmootherTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})

waLBerla_compile_test(FILES N1E1GMGhIndependenceTest.cpp common.hpp common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME N1E1GMGhIndependenceTest LABELS longrun LABELS n1e1)

waLBerla_compile_test(FILES N1E1GMGConvergenceTest.cpp common.hpp common.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME N1E1GMGConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1GMGConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

waLBerla_compile_test(FILES N1E1BasisTransformationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME N1E1BasisTransformationTest LABELS n1e1)

