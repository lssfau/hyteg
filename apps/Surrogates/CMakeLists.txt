waLBerla_add_executable(NAME compareP1Stencils
      FILES compareP1Stencils.cpp
      DEPENDS hyteg constant_stencil_operator)

waLBerla_add_executable( NAME PolynomialBlending
      FILES PolynomialBlending.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME VariableCoeffLaplace_Surrogates
      FILES VariableCoeffLaplace_Surrogates.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME Stokes2d_Surrogates
      FILES Stokes2d_Surrogates.cpp
      DEPENDS hyteg mixed_operator)

waLBerla_add_executable( NAME jumpCoefficient
      FILES jumpCoefficient.cpp
      DEPENDS hyteg mixed_operator)
