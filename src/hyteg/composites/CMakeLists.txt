target_sources( hyteg
  PRIVATE
  CCRStokesFunction.hpp
  CCRStokesOperator.hpp
  ConcatenatedOperator.hpp
  P0P1HelperFunctions.hpp
  P0P1UpwindOperator.hpp
  P1DGEP0StokesFunction.hpp
  P1DGEP0StokesOperator.hpp
  P1DGEP1StokesFunction.cpp
  P1P0StokesFunction.hpp
  P1P0StokesOperator.hpp
  P1P1UzawaDampingFactorEstimationOperator.hpp
  P1StokesFunction.hpp
  P2P1BlendingTaylorHoodStokesOperator.hpp
  P2P1SurrogateTaylorHoodStokesOperator.hpp
  P2P1TaylorHoodBlockFunction.hpp
  P2P1TaylorHoodFunction.hpp
  P2P2StabilizedStokesOperator.hpp
  P2P2StokesFunction.hpp
  StokesOperatorTraits.hpp
  StrongFreeSlipWrapper.hpp
  UnsteadyDiffusion.hpp
)

add_subdirectory( transport )
