waLBerla_link_files_to_builddir(*.msh)

waLBerla_compile_test(FILES P1InterpolateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME vCycleTest1 COMMAND $<TARGET_FILE:vCycleTest> )

waLBerla_compile_test(FILES Indexing/P1EdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/P1FaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleEdgeIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleFaceIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleFaceIndexingTest)

waLBerla_compile_test(FILES P1EnumerateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1EnumerateTest)

waLBerla_compile_test(FILES BubbleEnumerateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleEnumerateTest)

waLBerla_compile_test(FILES P1PackInfoTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1BPackInfoTest)

waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveStorageTest1 COMMAND $<TARGET_FILE:PrimitiveStorageTest> )
waLBerla_execute_test(NAME PrimitiveStorageTest3 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 3 )
waLBerla_execute_test(NAME PrimitiveStorageTest8 COMMAND $<TARGET_FILE:PrimitiveStorageTest> PROCESSES 8 )

waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BufferedCommunicationTest1 COMMAND $<TARGET_FILE:BufferedCommunicationTest> )
waLBerla_execute_test(NAME BufferedCommunicationTest3 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 3 )
waLBerla_execute_test(NAME BufferedCommunicationTest8 COMMAND $<TARGET_FILE:BufferedCommunicationTest> PROCESSES 8 )

waLBerla_compile_test(FILES P1DataHandlingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1DataHandlingTest)

waLBerla_compile_test(FILES PrimitiveMigrationTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveMigrationTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> )
waLBerla_execute_test(NAME PrimitiveMigrationTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 3 )
waLBerla_execute_test(NAME PrimitiveMigrationTest8 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 8 )

if( WALBERLA_BUILD_WITH_PARMETIS )
  waLBerla_compile_test(FILES ParMetisTestHHG.cpp DEPENDS tinyhhg_core core fmt)
  waLBerla_execute_test(NAME ParMetisTestHHG1 COMMAND $<TARGET_FILE:ParMetisTestHHG> )
  waLBerla_execute_test(NAME ParMetisTestHHG3 COMMAND $<TARGET_FILE:ParMetisTestHHG> PROCESSES 3 )
  waLBerla_execute_test(NAME ParMetisTestHHG8 COMMAND $<TARGET_FILE:ParMetisTestHHG> PROCESSES 8 )
endif()

waLBerla_compile_test(FILES MeshInfoTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME MeshInfoTest1 COMMAND $<TARGET_FILE:MeshInfoTest> )
waLBerla_execute_test(NAME MeshInfoTest2 COMMAND $<TARGET_FILE:MeshInfoTest> PROCESSES 2 )

waLBerla_compile_test(FILES FunctionSpaceDataTypesTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME FunctionSpaceDataTypesTest1 COMMAND $<TARGET_FILE:FunctionSpaceDataTypesTest> )