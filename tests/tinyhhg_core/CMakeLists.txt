waLBerla_link_files_to_builddir(*.msh)
waLBerla_link_files_to_builddir(../../data/param/*.prm)

waLBerla_compile_test(FILES p1integrateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME p1integrateTest)

waLBerla_compile_test(FILES fmgTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME fmgTest)

waLBerla_compile_test(FILES p1BubbleFaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME p1BubbleFaceIndexingTest)

waLBerla_compile_test(FILES p1BubbleEdgeIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME p1BubbleEdgeIndexingTest)

waLBerla_compile_test(FILES p1BubbleCommTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME p1BubbleCommTest)

waLBerla_compile_test(FILES P1BubblePackTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME P1BubblePackTest)

waLBerla_compile_test(FILES PrimitiveIDTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveIDTest)

waLBerla_compile_test(FILES PrimitiveDataTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveDataTest)

waLBerla_compile_test(FILES PrimitiveStorageTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME PrimitiveStorageTest)

waLBerla_compile_test(FILES BufferedCommunicationTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME BufferedCommunicationTest)