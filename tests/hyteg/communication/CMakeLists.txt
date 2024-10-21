waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BufferedCommunicationTest1 COMMAND $<TARGET_FILE:BufferedCommunicationTest>)
waLBerla_execute_test(NAME BufferedCommunicationTest3 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 3)
waLBerla_execute_test(NAME BufferedCommunicationTest8 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 8)

waLBerla_compile_test(FILES SyncingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SyncingTest1 COMMAND $<TARGET_FILE:SyncingTest>)
waLBerla_execute_test(NAME SyncingTest3 COMMAND $<TARGET_FILE:SyncingTest> PROCESSES 3)

waLBerla_compile_test(FILES CollectiveOperationsTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME CollectiveOperationsTest1 COMMAND $<TARGET_FILE:CollectiveOperationsTest>)
waLBerla_execute_test(NAME CollectiveOperationsTest3 COMMAND $<TARGET_FILE:CollectiveOperationsTest> PROCESSES 3)

