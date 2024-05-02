waLBerla_compile_test(FILES SphereToolsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SphereToolsTest1 COMMAND $<TARGET_FILE:SphereToolsTest>)
waLBerla_execute_test(NAME SphereToolsTest2 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 2)
waLBerla_execute_test(NAME SphereToolsTest4 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 4)

