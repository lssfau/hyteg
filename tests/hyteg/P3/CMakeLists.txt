waLBerla_add_test_executable( P3FunctionTest P3FunctionTest.cpp )
target_link_libraries       ( P3FunctionTest hyteg walberla::core )
waLBerla_execute_test(NAME P3FunctionTest1 COMMAND $<TARGET_FILE:P3FunctionTest> PROCESSES 1)
waLBerla_execute_test(NAME P3FunctionTest2 COMMAND $<TARGET_FILE:P3FunctionTest> PROCESSES 2)
