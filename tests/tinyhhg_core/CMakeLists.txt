waLBerla_link_files_to_builddir(*.msh)

waLBerla_compile_test(FILES p1integrateTest.cpp DEPENDS tinyhhg_core core fmt)
waLBerla_execute_test(NAME p1integrateTest)