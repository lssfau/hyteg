waLBerla_add_test_executable( FunctionWrapperTest FunctionWrapperTest.cpp )
target_link_libraries       ( FunctionWrapperTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionWrapperTest)
waLBerla_execute_test(NAME FunctionWrapperTestMPI COMMAND $<TARGET_FILE:FunctionWrapperTest> PROCESSES 2)

waLBerla_add_test_executable( BlockFunctionBasicTest BlockFunctionBasicTest.cpp )
target_link_libraries       ( BlockFunctionBasicTest hyteg walberla::core )
waLBerla_execute_test(NAME BlockFunctionBasicTest)

