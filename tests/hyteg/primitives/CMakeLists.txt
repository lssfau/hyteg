waLBerla_compile_test(FILES Face+CellGeometryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME Face+CellGeometryTest)

waLBerla_compile_test(FILES PrimitiveID+TypeTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveID+TypeTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES rootOnlySetupStorageTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME rootOnlySetupStorageTest COMMAND $<TARGET_FILE:rootOnlySetupStorageTest>)
waLBerla_execute_test(NAME rootOnlySetupStorageTest3 COMMAND $<TARGET_FILE:rootOnlySetupStorageTest> PROCESSES 3)
waLBerla_execute_test(NAME rootOnlySetupStorageTest8 COMMAND $<TARGET_FILE:rootOnlySetupStorageTest> PROCESSES 8)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageTest1 COMMAND $<TARGET_FILE:PrimitiveStorageTest>)
waLBerla_execute_test(NAME PrimitiveStorageTest3 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageTest8 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 8)

waLBerla_compile_test(FILES PrimitiveRelativeOrderingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveRelativeOrderingTest)

if (WALBERLA_BUILD_WITH_MPI)
    waLBerla_compile_test(FILES PrimitiveStorageParallelSetupTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest1 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest>)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest2 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 2)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest3 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 3)
    waLBerla_execute_test(NAME PrimitiveStorageParallelSetupTest4 COMMAND $<TARGET_FILE:PrimitiveStorageParallelSetupTest> PROCESSES 4)
endif()

