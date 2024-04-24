waLBerla_compile_test(FILES P2LinearCombinationFormTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2LinearCombinationFormTest)

waLBerla_compile_test(FILES RowSumTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME RowSumTest)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES HytegVsFenicsFormTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME HytegVsFenicsFormTest)
endif()

waLBerla_compile_test(FILES SingleRowIntegrationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME SingleRowIntegrationTest)

