target_sources( hyteg
    PRIVATE
    P2ElementwiseOperator.hpp
    P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp
    P1ToP2ElementwiseOperator.cpp
    P1ElementwiseOperator.hpp
    P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner.hpp
    DiagonalNonConstantOperator.hpp
    P2ToP1ElementwiseOperator.hpp
    P2P1ElementwiseBlendingStokesOperator.hpp
    P2ToP1ElementwiseOperator.cpp
    P1P1ElementwiseAffineEpsilonStokesOperator.hpp
    P2P1ElementwiseConstantCoefficientStokesOperator.hpp
    P2P1ElementwiseAffineEpsilonStokesOperator.hpp
    P2ElementwiseOperator.cpp
    P1ToP2ElementwiseOperator.hpp
    P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner.hpp
    P1ElementwiseOperator.cpp
    DiagonalNonConstantOperator.cpp     
)

