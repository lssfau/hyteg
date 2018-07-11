waLBerla_add_executable( NAME CompareInterpolate
        FILES CompareInterpolate.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES IndexFunctionVSManualVersion.cpp
        DEPENDS tinyhhg_core)

if( HHG_BUILD_WITH_PETSC )
    waLBerla_add_executable( NAME PetscMatVecMult
            FILES PetscMatVecMult.cpp
            DEPENDS tinyhhg_core core)
endif()

add_subdirectory(P1Benchmark)
add_subdirectory(P1CGBenchmark)