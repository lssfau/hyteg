waLBerla_compile_test(FILES P0FunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME P0FunctionTest1 COMMAND $<TARGET_FILE:P0FunctionTest>)
waLBerla_execute_test(NAME P0FunctionTest2 COMMAND $<TARGET_FILE:P0FunctionTest> PROCESSES 2)

waLBerla_compile_test(FILES P0AveragingTest.cpp DEPENDS hyteg core)
