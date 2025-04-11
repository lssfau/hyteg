add_library(mixed_operator)
target_link_libraries(mixed_operator PRIVATE hyteg core constant_stencil_operator)

target_sources(mixed_operator
      PRIVATE
      EGConvTestUtils.hpp
      EGOperators.hpp
      EGOperatorsNitscheBC.hpp
      HybridSmoother.cpp
      HybridSmoother.hpp
      P1BlendingStokesOperator.hpp
      P1P1ElementwiseAffineEpsilonStokesOperator.hpp
      P1P1StokesOperator.hpp
      P1StokesBlockPreconditioner.hpp
      P1ToP2ConstantOperator.hpp
      P2P1TaylorHoodStokesBlockPreconditioner.hpp
      P2P1TaylorHoodStokesOperator.hpp
      P2P1UzawaDampingFactorEstimationOperator.hpp
      P2P2UnstableStokesOperator.hpp
      P2ToP1ConstantOperator.hpp
      ScalarToVectorOperator.hpp
      VectorLaplaceOperator.cpp
      VectorLaplaceOperator.hpp
      VectorMassOperator.cpp
      VectorMassOperator.hpp
      VectorToScalarOperator.hpp
      VectorToVectorOperator.hpp
)
