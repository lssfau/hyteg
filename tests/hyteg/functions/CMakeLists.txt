waLBerla_add_test_executable( FunctionMemoryAllocationTest FunctionMemoryAllocationTest.cpp )
target_link_libraries       ( FunctionMemoryAllocationTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionMemoryAllocationTest1 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest>)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest2 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest8 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 8)

waLBerla_add_test_executable( FunctionSpaceDataTypesTest FunctionSpaceDataTypesTest.cpp )
target_link_libraries       ( FunctionSpaceDataTypesTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionSpaceDataTypesTest1 COMMAND $<TARGET_FILE:FunctionSpaceDataTypesTest>)

waLBerla_add_test_executable( FunctionIteratorTest FunctionIteratorTest.cpp )
target_link_libraries       ( FunctionIteratorTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionIteratorTest1 COMMAND $<TARGET_FILE:FunctionIteratorTest>)
waLBerla_execute_test(NAME FunctionIteratorTest2 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionIteratorTest8 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 8)

waLBerla_add_test_executable( FunctionMultiStoreTest FunctionMultiStoreTest.cpp )
target_link_libraries       ( FunctionMultiStoreTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionMultiStoreTest)

waLBerla_add_test_executable( FunctionRuleOfThreeTest FunctionRuleOfThreeTest.cpp )
target_link_libraries       ( FunctionRuleOfThreeTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionRuleOfThreeTest1 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest>)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest2 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest8 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 8)

waLBerla_add_test_executable( FEFunctionRegistryTest FEFunctionRegistryTest.cpp )
target_link_libraries       ( FEFunctionRegistryTest hyteg walberla::core )
waLBerla_execute_test(NAME FEFunctionRegistryTest COMMAND $<TARGET_FILE:FEFunctionRegistryTest>)

waLBerla_add_test_executable( FunctionMemorySerializationTest FunctionMemorySerializationTest.cpp )
target_link_libraries       ( FunctionMemorySerializationTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest>)
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8)

# Test some basic auxiliary stuff for application layer
waLBerla_add_test_executable( FindMaxMinMagTest FindMaxMinMagTest.cpp )
target_link_libraries       ( FindMaxMinMagTest hyteg walberla::core )
waLBerla_execute_test(NAME FindMaxMinMagTest1 COMMAND $<TARGET_FILE:FindMaxMinMagTest>)
waLBerla_execute_test(NAME FindMaxMinMagTest2 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 2)

waLBerla_add_test_executable( FunctionPropertiesTest FunctionPropertiesTest.cpp )
target_link_libraries       ( FunctionPropertiesTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionPropertiesTest)

waLBerla_add_test_executable( FunctionMultElementwiseTest FunctionMultElementwiseTest.cpp )
target_link_libraries       ( FunctionMultElementwiseTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionMultElementwiseTest)

waLBerla_add_test_executable( FunctionExtendedExpressionInterpolationTest FunctionExtendedExpressionInterpolationTest.cpp )
target_link_libraries       ( FunctionExtendedExpressionInterpolationTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionExtendedExpressionInterpolationTest)

waLBerla_add_test_executable( FunctionInterpolateOnceTest FunctionInterpolateOnceTest.cpp )
target_link_libraries       ( FunctionInterpolateOnceTest hyteg walberla::core )
waLBerla_execute_test(NAME FunctionInterpolateOnceTest1 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 1)
waLBerla_execute_test(NAME FunctionInterpolateOnceTest2 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionInterpolateOnceTest4 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 4)
