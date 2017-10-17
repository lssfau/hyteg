waLBerla_link_files_to_builddir(*.msh)


## Indexing ##

waLBerla_compile_test(FILES Indexing/P1EdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/P1FaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroEdgeIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME EdgeDoFOnMacroEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/EdgeDoFOnMacroFaceIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME EdgeDoFOnMacroFaceIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleEdgeIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleFaceIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleFaceIndexingTest)

waLBerla_compile_test(FILES Indexing/DGIndexingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME DGIndexingTest)

if(CMAKE_BUILD_TYPE MATCHES "Release")
  if(NOT WALBERLA_CXX_COMPILER_IS_MSVC )
    waLBerla_compile_test(FILES Indexing/ConstexprIndexTest.cpp DEPENDS tinyhhg_core core fmt)
    waLBerla_execute_test(NAME ConstexprIndexTest)
  endif()
endif()

###############


waLBerla_compile_test(FILES P1InterpolateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES DGInterpolateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME DGInterpolateTest)

waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME vCycleTest1 COMMAND $<TARGET_FILE:vCycleTest> )

waLBerla_compile_test(FILES P1CommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES BubbleCommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleCommTest)

waLBerla_compile_test(FILES DGCommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME DGCommTest)

waLBerla_compile_test(FILES DGEnumerateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME DGEnumerateTest)

waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES LevelinfoTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME LevelinfoTest)

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

if( HHG_BUILD_WITH_PETSC )
  waLBerla_compile_test(FILES PetscTest.cpp DEPENDS tinyhhg_core core fmt)
  waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest> )
  waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2 )
endif()

waLBerla_compile_test(FILES FunctionMemorySerializationTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME FunctionMemorySerializationTest1 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> )
waLBerla_execute_test(NAME FunctionMemorySerializationTest2 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 2 )
waLBerla_execute_test(NAME FunctionMemorySerializationTest8 COMMAND $<TARGET_FILE:FunctionMemorySerializationTest> PROCESSES 8 )

