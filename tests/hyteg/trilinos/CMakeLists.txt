if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES TrilinosVersionTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME TrilinosVersionTest)
endif ()

if (HYTEG_BUILD_WITH_TRILINOS)
    waLBerla_compile_test(FILES TrilinosSparseMatrixTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME TrilinosSparseMatrixTest)

    waLBerla_compile_test(FILES TrilinosDirectSolverTest.cpp DEPENDS hyteg core mixed_operator)
    #TODO: add one the solver is supported by mumps verison in image
    #waLBerla_execute_test(NAME TrilinosDirectSolverTest)

    waLBerla_compile_test(FILES TrilinosVectorTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME TrilinosVectorTest)
endif ()

