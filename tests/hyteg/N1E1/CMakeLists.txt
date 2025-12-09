# For some reason executing some tests takes forever with clang and gcc in Debug mode.
set(N1E1_CONDITIONAL_LONGRUN "")
if (     CMAKE_BUILD_TYPE MATCHES "Debug"
      AND (WALBERLA_CXX_COMPILER_IS_CLANG OR WALBERLA_CXX_COMPILER_IS_GNU)
)
    set(N1E1_CONDITIONAL_LONGRUN "longrun")
endif ()

waLBerla_add_test_executable( N1E1AddTest N1E1AddTest.cpp )
target_link_libraries       ( N1E1AddTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1AddTest LABELS n1e1)

waLBerla_add_test_executable( N1E1AssignTest N1E1AssignTest.cpp )
target_link_libraries       ( N1E1AssignTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1AssignTest LABELS n1e1)

waLBerla_add_test_executable( N1E1EvaluateTest N1E1EvaluateTest.cpp )
target_link_libraries       ( N1E1EvaluateTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1EvaluateTest LABELS n1e1)

waLBerla_add_test_executable( N1E1InterpolateTest N1E1InterpolateTest.cpp )
target_link_libraries       ( N1E1InterpolateTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1InterpolateTest LABELS n1e1)

waLBerla_add_test_executable( N1E1AdditiveCommunicationTest N1E1AdditiveCommunicationTest.cpp )
target_link_libraries       ( N1E1AdditiveCommunicationTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1AdditiveCommunicationTest LABELS n1e1)

waLBerla_add_test_executable( N1E1CurlCurlNullSpaceTest N1E1CurlCurlNullSpaceTest.cpp )
target_link_libraries       ( N1E1CurlCurlNullSpaceTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1CurlCurlNullSpaceTest LABELS n1e1)

waLBerla_add_test_executable( N1E1MassTest N1E1MassTest.cpp )
target_link_libraries       ( N1E1MassTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1MassTest LABELS n1e1)

waLBerla_add_test_executable( N1E1FormCurlCurlTest N1E1FormCurlCurlTest.cpp )
target_link_libraries       ( N1E1FormCurlCurlTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1FormCurlCurlTest LABELS n1e1)

waLBerla_add_test_executable( N1E1CurlCurlConvergenceTest N1E1CurlCurlConvergenceTest.cpp common.hpp common.cpp )
target_link_libraries       ( N1E1CurlCurlConvergenceTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1CurlCurlConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1CurlCurlConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

waLBerla_add_test_executable( N1E1ProlongationTest N1E1ProlongationTest.cpp )
target_link_libraries       ( N1E1ProlongationTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1ProlongationTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( N1E1RestrictionTest N1E1RestrictionTest.cpp )
    target_link_libraries       ( N1E1RestrictionTest hyteg walberla::core )
    waLBerla_execute_test(NAME N1E1RestrictionTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

waLBerla_add_test_executable( P1toN1E1GradientTest P1toN1E1GradientTest.cpp )
target_link_libraries       ( P1toN1E1GradientTest hyteg walberla::core )
waLBerla_execute_test(NAME P1toN1E1GradientTest LABELS n1e1)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( N1E1toP1LiftingTest N1E1toP1LiftingTest.cpp )
    target_link_libraries       ( N1E1toP1LiftingTest hyteg walberla::core )
    waLBerla_execute_test(NAME N1E1toP1LiftingTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})
endif ()

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( N1E1InverseDiagonalTest N1E1InverseDiagonalTest.cpp )
    target_link_libraries       ( N1E1InverseDiagonalTest hyteg walberla::core )
    waLBerla_execute_test(NAME N1E1InverseDiagonalTest LABELS n1e1)
endif ()

waLBerla_add_test_executable( N1E1BlendingTest N1E1BlendingTest.cpp )
target_link_libraries       ( N1E1BlendingTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1BlendingTest LABELS n1e1)

waLBerla_add_test_executable( N1E1LinearFormBlendingTest N1E1LinearFormBlendingTest.cpp )
target_link_libraries       ( N1E1LinearFormBlendingTest hyteg walberla::core )
waLBerla_execute_test(NAME N1E1LinearFormBlendingTest LABELS n1e1)

waLBerla_add_test_executable( HybridSmootherTest HybridSmootherTest.cpp common.hpp common.cpp )
target_link_libraries       ( HybridSmootherTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME HybridSmootherTest LABELS n1e1 ${N1E1_CONDITIONAL_LONGRUN})

waLBerla_add_test_executable( N1E1GMGhIndependenceTest N1E1GMGhIndependenceTest.cpp common.hpp common.cpp )
target_link_libraries       ( N1E1GMGhIndependenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME N1E1GMGhIndependenceTest LABELS longrun LABELS n1e1)

waLBerla_add_test_executable( N1E1GMGConvergenceTest N1E1GMGConvergenceTest.cpp common.hpp common.cpp )
target_link_libraries       ( N1E1GMGConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME N1E1GMGConvergenceTestLongrun COMMAND $<TARGET_FILE:N1E1GMGConvergenceTest> LABELS n1e1 longrun PROCESSES 4)

waLBerla_add_test_executable( N1E1BasisTransformationTest N1E1BasisTransformationTest.cpp )
target_link_libraries       ( N1E1BasisTransformationTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME N1E1BasisTransformationTest LABELS n1e1)

