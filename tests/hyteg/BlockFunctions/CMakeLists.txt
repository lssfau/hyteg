waLBerla_compile_test(FILES FunctionWrapperTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionWrapperTest)
waLBerla_execute_test(NAME FunctionWrapperTestMPI COMMAND $<TARGET_FILE:FunctionWrapperTest> PROCESSES 2)

waLBerla_compile_test(FILES BlockFunctionBasicTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BlockFunctionBasicTest)

