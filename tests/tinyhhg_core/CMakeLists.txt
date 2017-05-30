waLBerla_link_files_to_builddir(*.msh)
waLBerla_link_files_to_builddir(../../data/param/*.prm)

waLBerla_compile_test(FILES p1integrateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME p1integrateTest)

waLBerla_compile_test(FILES fmgTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME fmgTest)

waLBerla_compile_test(FILES p1BubbleFaceIndexingTest.cpp DEPENDS tinyhhg_core core)
waLBerla_execute_test(NAME p1BubbleFaceIndexingTest)