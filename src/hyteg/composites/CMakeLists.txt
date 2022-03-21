target_sources( hyteg
    PRIVATE
    P2P1TaylorHoodStokesBlockPreconditioner.hpp
    P1StokesBlockPreconditioner.hpp
    P2P1TaylorHoodBlockFunction.hpp
    P2P1TaylorHoodFunction.hpp
    P2P1SurrogateTaylorHoodStokesOperator.hpp
    P1P1StokesOperator.hpp
    P2P2UnstableStokesOperator.hpp
    P2P1TaylorHoodStokesOperator.hpp
    P2P1UzawaDampingFactorEstimationOperator.hpp
    P2P1BlendingTaylorHoodStokesOperator.hpp
    UnsteadyDiffusion.hpp
    P1Transport.hpp
    P2P2StokesFunction.hpp
    P1BlendingStokesOperator.hpp
    StokesOperatorTraits.hpp
    StrongFreeSlipWrapper.hpp
    ConcatenatedOperator.hpp
    P2P2StabilizedStokesOperator.hpp
    P1P1UzawaDampingFactorEstimationOperator.hpp
    P1EpsilonStokesOperator.hpp
    P1StokesFunction.hpp     
)
add_subdirectory( transport )
