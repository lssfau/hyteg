add_subdirectory(operators)

add_executable       ( compareP1Stencils compareP1Stencils.cpp )
target_link_libraries( compareP1Stencils hyteg constant_stencil_operator )

add_executable       ( PolynomialBlending PolynomialBlending.cpp )
target_link_libraries( PolynomialBlending hyteg mixed_operator )

add_executable       ( VariableCoeffLaplace_Surrogates VariableCoeffLaplace_Surrogates.cpp )
target_link_libraries( VariableCoeffLaplace_Surrogates hyteg mixed_operator )

add_executable       ( Stokes2d_Surrogates Stokes2d_Surrogates.cpp )
target_link_libraries( Stokes2d_Surrogates hyteg mixed_operator )

add_executable       ( jumpCoefficient jumpCoefficient.cpp )
target_link_libraries( jumpCoefficient hyteg mixed_operator )

if(HYTEG_BUILD_WITH_AVX AND WALBERLA_DOUBLE_ACCURACY)
      add_executable       ( surrogateBenchmark surrogateBenchmark.cpp )
      target_link_libraries( surrogateBenchmark hyteg opgen-mogli-diffusion constant_stencil_operator )
endif()
