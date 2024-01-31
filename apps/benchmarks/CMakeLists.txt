waLBerla_add_executable( NAME CompareInterpolate
        FILES CompareInterpolate.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES IndexFunctionVSManualVersion.cpp
        DEPENDS hyteg constantStencilOperator)

add_subdirectory(KernelBench)
if (WALBERLA_DOUBLE_ACCURACY)
    add_subdirectory(ApplyPerformanceAnalysis-2D-P2)
    add_subdirectory(ElementwiseOps)
    add_subdirectory(PetscCompare-2D-P2-Apply)
    add_subdirectory(P1Benchmark)
    add_subdirectory(P2Benchmark)
    add_subdirectory(P2OperatorBenchmarks)
    add_subdirectory(SnoopFilterIssueBenchmark)
endif ()

add_subdirectory(P2SolverBenchmark)
add_subdirectory(P1CGBenchmark)
add_subdirectory(ApplyBenchmark)

if( HYTEG_BUILD_WITH_PETSC )
    add_subdirectory(PetscCompare)
    add_subdirectory(PetscCompare-3D-P1-Apply)
endif()

add_subdirectory(surrogates)
