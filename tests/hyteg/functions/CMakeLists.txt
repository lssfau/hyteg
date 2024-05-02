waLBerla_compile_test(FILES FunctionMemoryAllocationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest1 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest>)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest2 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemoryAllocationTest8 COMMAND $<TARGET_FILE:FunctionMemoryAllocationTest> PROCESSES 8)

waLBerla_compile_test(FILES FunctionSpaceDataTypesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionSpaceDataTypesTest1 COMMAND $<TARGET_FILE:FunctionSpaceDataTypesTest>)

waLBerla_compile_test(FILES FunctionIteratorTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionIteratorTest1 COMMAND $<TARGET_FILE:FunctionIteratorTest>)
waLBerla_execute_test(NAME FunctionIteratorTest2 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionIteratorTest8 COMMAND $<TARGET_FILE:FunctionIteratorTest> PROCESSES 8)

waLBerla_compile_test(FILES FunctionMultiStoreTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMultiStoreTest)

waLBerla_compile_test(FILES FunctionRuleOfThreeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest1 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest>)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest2 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionRuleOfThreeTest8 COMMAND $<TARGET_FILE:FunctionRuleOfThreeTest> PROCESSES 8)

waLBerla_compile_test(FILES FEFunctionRegistryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FEFunctionRegistryTest COMMAND $<TARGET_FILE:FEFunctionRegistryTest>)

waLBerla_compile_test(FILES FunctionMemorySerializationTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest>)
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2)
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8)

# Test some basic auxiliary stuff for application layer
waLBerla_compile_test(FILES FindMaxMinMagTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FindMaxMinMagTest1 COMMAND $<TARGET_FILE:FindMaxMinMagTest>)
waLBerla_execute_test(NAME FindMaxMinMagTest2 COMMAND $<TARGET_FILE:FindMaxMinMagTest> PROCESSES 2)

waLBerla_compile_test(FILES FunctionPropertiesTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionPropertiesTest)

waLBerla_compile_test(FILES FunctionMultElementwiseTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionMultElementwiseTest)

waLBerla_compile_test(FILES FunctionExtendedExpressionInterpolationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME FunctionExtendedExpressionInterpolationTest)

