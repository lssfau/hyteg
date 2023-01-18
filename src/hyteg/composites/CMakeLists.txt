target_sources( hyteg
    PRIVATE
    ConcatenatedOperator.hpp
    P0P1UpwindOperator.hpp
    P1BlendingStokesOperator.hpp
    P1EpsilonStokesOperator.hpp
    P1P0StokesFunction.hpp
    P1P0StokesOperator.hpp
    P1P1StokesOperator.hpp
    P1P1UzawaDampingFactorEstimationOperator.hpp
    P1StokesBlockPreconditioner.hpp
    P1StokesFunction.hpp
    P1Transport.hpp
    P2P1BlendingTaylorHoodStokesOperator.hpp
    P2P1SurrogateTaylorHoodStokesOperator.hpp
    P2P1TaylorHoodBlockFunction.hpp
    P2P1TaylorHoodFunction.hpp
    P2P1TaylorHoodStokesBlockPreconditioner.hpp
    P2P1TaylorHoodStokesOperator.hpp
    P2P1UzawaDampingFactorEstimationOperator.hpp
    P2P2StabilizedStokesOperator.hpp
    P2P2StokesFunction.hpp
    P2P2UnstableStokesOperator.hpp
    StokesOperatorTraits.hpp
    StrongFreeSlipWrapper.hpp
    UnsteadyDiffusion.hpp
)
add_subdirectory( transport )
