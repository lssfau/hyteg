waLBerla_compile_test(FILES invertElementwiseTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME invertElementwiseTest)
waLBerla_execute_test(NAME invertElementwiseTestMPI COMMAND $<TARGET_FILE:invertElementwiseTest> PROCESSES 3)

waLBerla_compile_test(FILES CircularMapFPETest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME CircularMapFPETest)
waLBerla_execute_test(NAME CircularMapFPETestMPI COMMAND $<TARGET_FILE:CircularMapFPETest> PROCESSES 2)

