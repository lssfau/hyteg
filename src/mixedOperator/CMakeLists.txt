add_library(mixedOperator)
target_link_libraries(mixedOperator PRIVATE hyteg core constantStencilOperator)

target_sources(mixedOperator
      PRIVATE
      EGConvTestUtils.hpp
      EGOperators.hpp
      EGOperatorsNitscheBC.hpp
      HybridSmoother.hpp
      HybridSmoother.cpp
      P1BlendingStokesOperator.hpp
      P1P1ElementwiseAffineEpsilonStokesOperator.hpp
      P1ToP2ConstantOperator.hpp
      P2P1TaylorHoodStokesBlockPreconditioner.hpp
      P2P1TaylorHoodStokesOperator.hpp
      P2P1UzawaDampingFactorEstimationOperator.hpp
      P1P1StokesOperator.hpp
      P1StokesBlockPreconditioner.hpp
      P2ToP1ConstantOperator.hpp
      P2P2UnstableStokesOperator.hpp
      ScalarToVectorOperator.hpp
      VectorLaplaceOperator.cpp
      VectorLaplaceOperator.hpp
      VectorMassOperator.cpp
      VectorMassOperator.hpp
      VectorToScalarOperator.hpp
      VectorToVectorOperator.hpp
)