waLBerla_compile_test(FILES EdgeDoFToVertexDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFToVertexDoFOperatorTest)

waLBerla_compile_test(FILES GEMVSmokeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME GEMVSmokeTest)

waLBerla_compile_test(FILES SimpleVertexDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SimpleVertexDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES P1SurrogateOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1SurrogateOperatorTest)

waLBerla_compile_test(FILES P1ElementwiseSurrogateOperatorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P1ElementwiseSurrogateOperatorTest)
waLBerla_execute_test(NAME P1ElementwiseSurrogateOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 3)

waLBerla_compile_test(FILES DivergenceOperatorTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-divergence)
waLBerla_execute_test(NAME DivergenceOperatorTest)

waLBerla_compile_test(FILES ElementwiseOperatorAdditiveApplyTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME ElementwiseOperatorAdditiveApplyTest)

waLBerla_compile_test(FILES ElementwiseConstantAndDofValueGEMVMangagerTest.cpp DEPENDS hyteg core constant_stencil_operator mixed_operator elementwise_dof_value_operator)
waLBerla_execute_test(NAME ElementwiseConstantAndDofValueGEMVMangagerTest)

waLBerla_compile_test(FILES ElementwiseDoFValueOperatorTest.cpp DEPENDS hyteg core elementwise_dof_value_operator)
waLBerla_execute_test(NAME ElementwiseDoFValueOperatorTest)

waLBerla_compile_test(FILES ViscousOperatorsTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-viscousblock)
waLBerla_execute_test(NAME ViscousOperatorsTest)

waLBerla_compile_test(FILES VectorOperatorBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorOperatorBasicTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES ElementwiseOperatorPetscTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest>)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest2 COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest> PROCESSES 2)
endif ()

waLBerla_compile_test(FILES OperatorWrapperBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME OperatorWrapperBasicTest)

waLBerla_compile_test(FILES BlockOperatorBasicTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME BlockOperatorBasicTest)

waLBerla_compile_test(FILES BlockOperatorGSSmoothableTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME BlockOperatorGSSmoothableTest)

waLBerla_compile_test(FILES EdgeDoFToEdgeDoFOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME EdgeDoFToEdgeDoFOperatorTest)

waLBerla_compile_test(FILES MassOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME MassOperatorTest)
    waLBerla_execute_test(NAME MassOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 2)
endif ()

waLBerla_compile_test(FILES DiagonalNonConstantOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME DiagonalNonConstantOperatorTest)

waLBerla_compile_test(FILES P1ElementwiseOperatorTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1ElementwiseOperatorTest)

