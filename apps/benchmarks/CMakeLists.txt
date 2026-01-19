add_executable       ( CompareInterpolate CompareInterpolate.cpp )
target_link_libraries( CompareInterpolate hyteg )

add_executable       ( IndexFunctionVSManualVersion IndexFunctionVSManualVersion.cpp )
target_link_libraries( IndexFunctionVSManualVersion hyteg constant_stencil_operator )

add_executable       ( CollectiveOperationsBenchmark CollectiveOperationsBenchmark.cpp )
target_link_libraries( CollectiveOperationsBenchmark hyteg )

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
