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

    generated/FullStokesRotationWithFSPenalty/P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator.hpp
    generated/FullStokesRotationWithFSPenalty/P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator.cpp

    generated/Ones/P2Ones.hpp
    generated/Ones/P2Ones.cpp
)
