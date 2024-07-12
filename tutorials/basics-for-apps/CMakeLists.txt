waLBerla_add_executable( NAME BA.01_PrimitiveStorage
        FILES BA.01_PrimitiveStorage.cpp
        DEPENDS hyteg )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.01_PrimitiveStorage )

waLBerla_add_executable( NAME BA.02_FEFunctions
        FILES BA.02_FEFunctions.cpp
        DEPENDS hyteg )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.02_FEFunctions )

waLBerla_add_executable( NAME BA.03_ExportForVisualisation
        FILES BA.03_ExportForVisualisation.cpp
        DEPENDS hyteg )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.03_ExportForVisualisation )

add_subdirectory( BA.04_BoundaryConditions )

waLBerla_add_executable( NAME BA.05_AMR
        FILES BA.05_AMR.cpp
        DEPENDS hyteg )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.05_AMR )

waLBerla_add_executable( NAME BA.06_ErrorEstimator
        FILES BA.06_ErrorEstimator.cpp
        DEPENDS core hyteg opgen-mass constant_stencil_operator )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.06_ErrorEstimator )
