waLBerla_compile_test(FILES EigenBasicTest.cpp)
waLBerla_execute_test(NAME EigenBasicTest)
waLBerla_execute_test(NAME EigenBasicTestMPI COMMAND $<TARGET_FILE:EigenBasicTest> PROCESSES 2)

waLBerla_compile_test(FILES EigenAssemblyTest.cpp)
waLBerla_execute_test(NAME EigenAssemblyTest)

waLBerla_compile_test(FILES EigenSparseSolverTest.cpp)
waLBerla_execute_test(NAME EigenSparseSolverTest)

