waLBerla_add_test_executable( EdgeDoFToVertexDoFOperatorTest EdgeDoFToVertexDoFOperatorTest.cpp )
target_link_libraries       ( EdgeDoFToVertexDoFOperatorTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFToVertexDoFOperatorTest)

waLBerla_add_test_executable( GEMVSmokeTest GEMVSmokeTest.cpp )
target_link_libraries       ( GEMVSmokeTest hyteg walberla::core )
waLBerla_execute_test(NAME GEMVSmokeTest)

waLBerla_add_test_executable( SimpleVertexDoFToEdgeDoFOperatorTest SimpleVertexDoFToEdgeDoFOperatorTest.cpp )
target_link_libraries       ( SimpleVertexDoFToEdgeDoFOperatorTest hyteg walberla::core )
waLBerla_execute_test(NAME SimpleVertexDoFToEdgeDoFOperatorTest)

waLBerla_add_test_executable( P1SurrogateOperatorTest P1SurrogateOperatorTest.cpp )
target_link_libraries       ( P1SurrogateOperatorTest hyteg walberla::core )
waLBerla_execute_test(NAME P1SurrogateOperatorTest)

waLBerla_add_test_executable( P1ElementwiseSurrogateOperatorTest P1ElementwiseSurrogateOperatorTest.cpp )
target_link_libraries       ( P1ElementwiseSurrogateOperatorTest hyteg walberla::core )
waLBerla_execute_test(NAME P1ElementwiseSurrogateOperatorTest)
waLBerla_execute_test(NAME P1ElementwiseSurrogateOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 3)

waLBerla_add_test_executable( DivergenceOperatorTest DivergenceOperatorTest.cpp )
target_link_libraries       ( DivergenceOperatorTest hyteg walberla::core mixed_operator opgen-composites-divergence )
waLBerla_execute_test(NAME DivergenceOperatorTest)

waLBerla_add_test_executable( ElementwiseOperatorAdditiveApplyTest ElementwiseOperatorAdditiveApplyTest.cpp )
target_link_libraries       ( ElementwiseOperatorAdditiveApplyTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME ElementwiseOperatorAdditiveApplyTest)

waLBerla_add_test_executable( ElementwiseConstantAndDofValueGEMVMangagerTest ElementwiseConstantAndDofValueGEMVMangagerTest.cpp )
target_link_libraries       ( ElementwiseConstantAndDofValueGEMVMangagerTest hyteg walberla::core constant_stencil_operator mixed_operator elementwise_dof_value_operator )
waLBerla_execute_test(NAME ElementwiseConstantAndDofValueGEMVMangagerTest)

waLBerla_add_test_executable( ElementwiseDoFValueOperatorTest ElementwiseDoFValueOperatorTest.cpp )
target_link_libraries       ( ElementwiseDoFValueOperatorTest hyteg walberla::core elementwise_dof_value_operator )
waLBerla_execute_test(NAME ElementwiseDoFValueOperatorTest)

waLBerla_add_test_executable( ViscousOperatorsTest ViscousOperatorsTest.cpp )
target_link_libraries       ( ViscousOperatorsTest hyteg walberla::core mixed_operator opgen-composites-viscousblock )
waLBerla_execute_test(NAME ViscousOperatorsTest)

waLBerla_add_test_executable( VectorOperatorBasicTest VectorOperatorBasicTest.cpp )
target_link_libraries       ( VectorOperatorBasicTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME VectorOperatorBasicTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( ElementwiseOperatorPetscTest ElementwiseOperatorPetscTest.cpp )
    target_link_libraries       ( ElementwiseOperatorPetscTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest>)
    waLBerla_execute_test(NAME ElementwiseOperatorPetscTest2 COMMAND $<TARGET_FILE:ElementwiseOperatorPetscTest> PROCESSES 2)
endif ()

waLBerla_add_test_executable( OperatorWrapperBasicTest OperatorWrapperBasicTest.cpp )
target_link_libraries       ( OperatorWrapperBasicTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME OperatorWrapperBasicTest)

waLBerla_add_test_executable( BlockOperatorBasicTest BlockOperatorBasicTest.cpp )
target_link_libraries       ( BlockOperatorBasicTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME BlockOperatorBasicTest)

waLBerla_add_test_executable( BlockOperatorGSSmoothableTest BlockOperatorGSSmoothableTest.cpp )
target_link_libraries       ( BlockOperatorGSSmoothableTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME BlockOperatorGSSmoothableTest)

waLBerla_add_test_executable( EdgeDoFToEdgeDoFOperatorTest EdgeDoFToEdgeDoFOperatorTest.cpp )
target_link_libraries       ( EdgeDoFToEdgeDoFOperatorTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME EdgeDoFToEdgeDoFOperatorTest)

waLBerla_add_test_executable( MassOperatorTest MassOperatorTest.cpp )
target_link_libraries       ( MassOperatorTest hyteg walberla::core constant_stencil_operator )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME MassOperatorTest)
    waLBerla_execute_test(NAME MassOperatorTestMPI COMMAND $<TARGET_FILE:MassOperatorTest> PROCESSES 2)
endif ()

waLBerla_add_test_executable( DiagonalNonConstantOperatorTest DiagonalNonConstantOperatorTest.cpp )
target_link_libraries       ( DiagonalNonConstantOperatorTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME DiagonalNonConstantOperatorTest)

