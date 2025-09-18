target_sources( hyteg
    PRIVATE
    DivKGrad/P1ElementwiseDivKGrad_AnnulusMap_float64.hpp
    DivKGrad/P1ElementwiseDivKGrad_AnnulusMap_float64.cpp
    DivKGrad/P1ElementwiseDivKGrad_IcosahedralShellMap_float64.hpp
    DivKGrad/P1ElementwiseDivKGrad_IcosahedralShellMap_float64.cpp

    EpsilonRotation/P2VectorElementwiseEpsilonRotation_IcosahedralShellMap_fused_quadloops_float64.hpp
    EpsilonRotation/P2VectorElementwiseEpsilonRotation_IcosahedralShellMap_fused_quadloops_float64.cpp

    EpsilonRotation/P2VectorElementwiseEpsilonRotationAnnulusMap.hpp
    EpsilonRotation/P2VectorElementwiseEpsilonRotationAnnulusMap.cpp

    EpsilonRotation/P2VectorElementwiseEpsilonRotationIcosahedralShellMapVectorizedOperator.hpp
    EpsilonRotation/P2VectorElementwiseEpsilonRotationIcosahedralShellMapVectorizedOperator.cpp

    generated/FullStokesViscoplastic/P2VectorElementwiseFullStokesViscoplastic.hpp
    generated/FullStokesViscoplastic/P2VectorElementwiseFullStokesViscoplastic.cpp

    generated/EvaluateViscosityViscoplastic/P2EvaluateViscosityViscoplastic.hpp
    generated/EvaluateViscosityViscoplastic/P2EvaluateViscosityViscoplastic.cpp

    generated/RadialGradient/P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator.hpp
    generated/RadialGradient/P2ElementwiseRadialGradientSurfaceIcosahedralShellMapOperator.cpp

    generated/RadialGradient/P2ElementwiseRadialGradientSurfaceOperator.hpp
    generated/RadialGradient/P2ElementwiseRadialGradientSurfaceOperator.cpp

    generated/FullStokesRotationWithFSPenalty/P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator.hpp
    generated/FullStokesRotationWithFSPenalty/P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator.cpp

    generated/Ones/P2Ones.hpp
    generated/Ones/P2Ones.cpp

    generated/GradientRotation/P1ToP2VectorGradientRotationIcosahedralShellMapOperator.hpp
    generated/GradientRotation/P1ToP2VectorGradientRotationIcosahedralShellMapOperator.cpp

    generated/DivergenceRotation/P2VectorToP1DivergenceRotationIcosahedralShellMapOperator.hpp
    generated/DivergenceRotation/P2VectorToP1DivergenceRotationIcosahedralShellMapOperator.cpp

    generated/BoundaryMass/P2ElementwiseBoundaryMass.hpp
    generated/BoundaryMass/P2ElementwiseBoundaryMass.cpp

    generated/GradientBoundaryMass/P2ElementwiseGradientBoundaryMass.hpp
    generated/GradientBoundaryMass/P2ElementwiseGradientBoundaryMass.cpp

    generated/BoundaryMass/P2ElementwiseBoundaryMassIcosahedralShellMapOperator.hpp
    generated/BoundaryMass/P2ElementwiseBoundaryMassIcosahedralShellMapOperator.cpp

    generated/BoundaryKMass/P2ElementwiseBoundaryKMassIcosahedralShellMap.hpp
    generated/BoundaryKMass/P2ElementwiseBoundaryKMassIcosahedralShellMap.cpp

    generated/BoundaryDeltaMass/P2ElementwiseBoundaryDeltaMassIcosahedralShellMap.hpp
    generated/BoundaryDeltaMass/P2ElementwiseBoundaryDeltaMassIcosahedralShellMap.cpp

    generated/EpsilonRotation/P2VectorEpsilonRotationIcosahedralShellMapOperator.hpp
    generated/EpsilonRotation/P2VectorEpsilonRotationIcosahedralShellMapOperator.cpp

    generated/EpsilonRotationWithFSPenalty/P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator.hpp
    generated/EpsilonRotationWithFSPenalty/P2VectorElementwiseEpsilonRotationWithFSPenaltyIcosahedralShellMapOperator.cpp
)
