waLBerla_add_test_executable( IsoparametricPoissonConvergenceTest IsoparametricPoissonConvergenceTest.cpp )
target_link_libraries       ( IsoparametricPoissonConvergenceTest hyteg walberla::core mixed_operator opgen-diffusion opgen-mass )
waLBerla_execute_test(NAME IsoparametricPoissonConvergenceTest)