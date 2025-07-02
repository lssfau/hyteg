target_sources( hyteg
    PRIVATE
    StokesPCGSolver.hpp
    EmptySolver.hpp
    FlexibleMultigridSolver.hpp
    FullMultigridSolver.hpp
    Smoothables.hpp
    GeometricMultigridSolver.hpp
    SymmetricSORSmoother.hpp
    UzawaSmoother.hpp
    Solver.hpp
    GaussSeidelSmoother.hpp
    SymmetricGaussSeidelSmoother.hpp
    ChebyshevSmoother.hpp
    MinresSolver.hpp
    StokesPCGSolverOld.hpp
    FAS.hpp
    WeightedJacobiSmoother.hpp
    CGSolver.hpp
    SORSmoother.hpp     
    SubstitutePreconditioner.hpp
    ApplyInverseDiagonalWrapper.hpp
    GMRESSolver.hpp
    FGMRESSolver.hpp
    FunctionMultiplicationPreconditioner.hpp
)
add_subdirectory( solvertemplates )
add_subdirectory( preconditioners )
add_subdirectory( controlflow )
