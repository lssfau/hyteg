waLBerla_add_executable( NAME BA.01_PrimitiveStorage
        FILES BA.01_PrimitiveStorage.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME BA.05_AMR
        FILES BA.05_AMR.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME BA.06_ErrorEstimator
        FILES BA.06_ErrorEstimator.cpp
        DEPENDS core hyteg opgen-mass constant_stencil_operator)

add_subdirectory(BA.04_BoundaryConditions)
