waLBerla_add_test_executable( P0FunctionTest P0FunctionTest.cpp )
target_link_libraries       ( P0FunctionTest hyteg walberla::core )
waLBerla_execute_test(NAME P0FunctionTest1 COMMAND $<TARGET_FILE:P0FunctionTest>)
waLBerla_execute_test(NAME P0FunctionTest2 COMMAND $<TARGET_FILE:P0FunctionTest> PROCESSES 2)

waLBerla_add_test_executable( P0AveragingTest P0AveragingTest.cpp )
target_link_libraries       ( P0AveragingTest hyteg walberla::core )
waLBerla_execute_test(NAME P0AveragingTest1 COMMAND $<TARGET_FILE:P0AveragingTest>)
waLBerla_execute_test(NAME P0AveragingTest2 COMMAND $<TARGET_FILE:P0AveragingTest> PROCESSES 2)
