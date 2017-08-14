waLBerla_link_files_to_builddir(*.msh)
waLBerla_link_files_to_builddir(../../data/param/*.prm)

waLBerla_compile_test(FILES P1InterpolateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1InterpolateTest)

waLBerla_compile_test(FILES vCycleTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME vCycleTest)

waLBerla_compile_test(FILES Indexing/P1EdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1EdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/P1FaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME P1FaceIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleEdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME BubbleEdgeIndexingTest)

waLBerla_compile_test(FILES Indexing/BubbleFaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME BubbleFaceIndexingTest)

waLBerla_compile_test(FILES p1BubbleCommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME p1BubbleCommTest)

waLBerla_compile_test(FILES P1CommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1CommTest)

waLBerla_compile_test(FILES BubbleEnumerateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BubbleEnumerateTest)

waLBerla_compile_test(FILES P1PackInfoTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1BPackInfoTest)

waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME PrimitiveStorageTest)

waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME BufferedCommunicationTest)

waLBerla_compile_test(FILES p1BubbleVertexCommTest.cpp DEPENDS tinyhhg_core core fmt)

waLBerla_compile_test(FILES P1DataHandlingTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1DataHandlingTest)
