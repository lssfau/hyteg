waLBerla_add_test_executable( EigenBasicTest EigenBasicTest.cpp )
target_link_libraries       ( EigenBasicTest hyteg walberla::core )
waLBerla_execute_test(NAME EigenBasicTest)
waLBerla_execute_test(NAME EigenBasicTestMPI COMMAND $<TARGET_FILE:EigenBasicTest> PROCESSES 2)

waLBerla_add_test_executable( EigenAssemblyTest EigenAssemblyTest.cpp )
target_link_libraries       ( EigenAssemblyTest hyteg walberla::core )
waLBerla_execute_test(NAME EigenAssemblyTest)

waLBerla_add_test_executable( EigenSparseSolverTest EigenSparseSolverTest.cpp )
target_link_libraries       ( EigenSparseSolverTest hyteg walberla::core )
waLBerla_execute_test(NAME EigenSparseSolverTest)
