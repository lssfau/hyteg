target_sources( hyteg
    PRIVATE
    SchurCGCoarseGridSolver.hpp
    SchurCGOuterLoopSolver.hpp
    SchurCGPoissonNeumannSolver.hpp
    SchurChebyshevSmoother.hpp
    SchurFullMultigridSolver.hpp
    SchurIdentityPreconditioner.hpp
    SchurSolver.hpp
    SchurInverseDiagMassSolver.hpp
    SchurInverseLumpedMassSolver.hpp
    SchurMultigridSolver.hpp
    SchurCGMassSolver.hpp
    SchurBFBTSolver.hpp
    SchurPoissonNeumannChebyshevSmoother.hpp
    SchurPoissonNeumannMultigridSolver.hpp
)

