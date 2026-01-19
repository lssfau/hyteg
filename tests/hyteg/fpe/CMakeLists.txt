waLBerla_add_test_executable( invertElementwiseTest invertElementwiseTest.cpp )
target_link_libraries       ( invertElementwiseTest hyteg walberla::core )
waLBerla_execute_test(NAME invertElementwiseTest)
waLBerla_execute_test(NAME invertElementwiseTestMPI COMMAND $<TARGET_FILE:invertElementwiseTest> PROCESSES 3)

waLBerla_add_test_executable( CircularMapFPETest CircularMapFPETest.cpp )
target_link_libraries       ( CircularMapFPETest hyteg walberla::core )
waLBerla_execute_test(NAME CircularMapFPETest)
waLBerla_execute_test(NAME CircularMapFPETestMPI COMMAND $<TARGET_FILE:CircularMapFPETest> PROCESSES 2)

