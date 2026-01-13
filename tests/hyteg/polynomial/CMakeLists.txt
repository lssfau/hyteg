waLBerla_add_test_executable( PolynomialTest PolynomialTest.cpp )
target_link_libraries       ( PolynomialTest hyteg walberla::core )
waLBerla_execute_test(NAME PolynomialTest)

waLBerla_add_test_executable( StencilPolynomialTest StencilPolynomialTest.cpp )
target_link_libraries       ( StencilPolynomialTest hyteg walberla::core )
waLBerla_execute_test(NAME StencilPolynomialTest)

waLBerla_add_test_executable( PiecewiseLSQPInterpolationTest PiecewiseLSQPInterpolationTest.cpp )
target_link_libraries       ( PiecewiseLSQPInterpolationTest hyteg walberla::core )
waLBerla_execute_test(NAME PiecewiseLSQPInterpolationTest)

waLBerla_add_test_executable( LeastSquaresTest LeastSquaresTest.cpp )
target_link_libraries       ( LeastSquaresTest hyteg walberla::core )
waLBerla_execute_test(NAME LeastSquaresTest)
