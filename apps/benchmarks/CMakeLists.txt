waLBerla_add_executable( NAME CompareInterpolate
        FILES CompareInterpolate.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES IndexFunctionVSManualVersion.cpp
        DEPENDS tinyhhg_core)

add_subdirectory(P1Benchmark)
#add_subdirectory(P1CGBenchmark)

if( HHG_BUILD_WITH_PETSC )
    add_subdirectory(PetscCompare)
    add_subdirectory(PetscCompare-2D-P2-Apply)
    add_subdirectory(ApplyPerformanceAnalysis-2D-P2)
endif()