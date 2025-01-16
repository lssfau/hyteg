waLBerla_compile_test(FILES PolynomialTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PolynomialTest)

waLBerla_compile_test(FILES LeastSquaresTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME LeastSquaresTest)