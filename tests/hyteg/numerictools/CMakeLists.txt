waLBerla_add_test_executable( L2Test L2Test.cpp )
target_link_libraries       ( L2Test hyteg walberla::core )
waLBerla_execute_test(NAME L2Test)
waLBerla_execute_test(NAME L2TestMPI COMMAND $<TARGET_FILE:L2Test> PROCESSES 3)

waLBerla_add_test_executable( SpectrumEstimationTest SpectrumEstimationTest.cpp )
target_link_libraries       ( SpectrumEstimationTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME SpectrumEstimationTest)
waLBerla_execute_test(NAME SpectrumEstimationTestMPI COMMAND $<TARGET_FILE:SpectrumEstimationTest> PROCESSES 3)

