if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_add_test_executable( TrilinosVersionTest TrilinosVersionTest.cpp )
    target_link_libraries       ( TrilinosVersionTest hyteg walberla::core )
    waLBerla_execute_test(NAME TrilinosVersionTest)
endif ()

if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_add_test_executable( TrilinosSparseMatrixTest TrilinosSparseMatrixTest.cpp )
    target_link_libraries       ( TrilinosSparseMatrixTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME TrilinosSparseMatrixTest)

    waLBerla_add_test_executable( TrilinosDirectSolverTest TrilinosDirectSolverTest.cpp )
    target_link_libraries       ( TrilinosDirectSolverTest hyteg walberla::core mixed_operator )
    #TODO: add one the solver is supported by mumps verison in image
    #waLBerla_execute_test(NAME TrilinosDirectSolverTest)

    waLBerla_add_test_executable( TrilinosVectorTest TrilinosVectorTest.cpp )
    target_link_libraries       ( TrilinosVectorTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME TrilinosVectorTest)
endif ()

