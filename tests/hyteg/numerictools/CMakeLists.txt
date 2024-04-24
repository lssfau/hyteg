waLBerla_compile_test(FILES L2Test.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME L2Test)
waLBerla_execute_test(NAME L2TestMPI COMMAND $<TARGET_FILE:L2Test> PROCESSES 3)

waLBerla_compile_test(FILES SpectrumEstimationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME SpectrumEstimationTest)
waLBerla_execute_test(NAME SpectrumEstimationTestMPI COMMAND $<TARGET_FILE:SpectrumEstimationTest> PROCESSES 3)

