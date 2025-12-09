waLBerla_add_test_executable( MatrixConstructFromEigenTest MatrixConstructFromEigenTest.cpp )
target_link_libraries       ( MatrixConstructFromEigenTest hyteg walberla::core )
waLBerla_execute_test(NAME MatrixConstructFromEigenTest)

