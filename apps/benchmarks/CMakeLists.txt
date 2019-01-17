waLBerla_add_executable( NAME CompareInterpolate
        FILES CompareInterpolate.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES IndexFunctionVSManualVersion.cpp
        DEPENDS tinyhhg_core)

add_subdirectory(P1Benchmark)
add_subdirectory(P2Benchmark)
add_subdirectory(P1CGBenchmark)
add_subdirectory(ApplyPerformanceAnalysis-2D-P2)
add_subdirectory(ApplyBenchmark)

if( HHG_BUILD_WITH_PETSC )
    add_subdirectory(PetscCompare)
    add_subdirectory(PetscCompare-2D-P2-Apply)
    add_subdirectory(PetscCompare-3D-P1-Apply)
endif()