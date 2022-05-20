waLBerla_add_executable( NAME CompareInterpolate
        FILES CompareInterpolate.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES IndexFunctionVSManualVersion.cpp
        DEPENDS hyteg)

add_subdirectory(ElementwiseOps)
add_subdirectory(KernelBench)
add_subdirectory(P1Benchmark)
add_subdirectory(P2Benchmark)
add_subdirectory(P2OperatorBenchmarks)
add_subdirectory(P2SolverBenchmark)
add_subdirectory(P1CGBenchmark)
add_subdirectory(ApplyPerformanceAnalysis-2D-P2)
add_subdirectory(ApplyBenchmark)
add_subdirectory(SnoopFilterIssueBenchmark)
add_subdirectory(PetscCompare-2D-P2-Apply)

if( HYTEG_BUILD_WITH_PETSC )
    add_subdirectory(PetscCompare)
    add_subdirectory(PetscCompare-3D-P1-Apply)
endif()

add_subdirectory(surrogates)
