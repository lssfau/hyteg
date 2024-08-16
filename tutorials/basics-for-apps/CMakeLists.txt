waLBerla_add_executable( NAME BA.01_PrimitiveStorage
        FILES BA.01_PrimitiveStorage.cpp
        DEPENDS hyteg )
if ( CMAKE_BUILD_TYPE MATCHES "Release" )
   waLBerla_execute_test( NO_MODULE_LABEL NAME BA.01_PrimitiveStorage )
endif ()

waLBerla_add_executable( NAME BA.02_FEFunctions
        FILES BA.02_FEFunctions.cpp
        DEPENDS hyteg )
if( CMAKE_BUILD_TYPE MATCHES "Release" )
   waLBerla_execute_test( NO_MODULE_LABEL NAME BA.02_FEFunctions )
endif()

waLBerla_add_executable( NAME BA.03_ExportForVisualisation
        FILES BA.03_ExportForVisualisation.cpp
        DEPENDS hyteg )
if( CMAKE_BUILD_TYPE MATCHES "Release" )
   waLBerla_execute_test( NO_MODULE_LABEL NAME BA.03_ExportForVisualisation )
endif()

add_subdirectory( BA.04_BoundaryConditions )

waLBerla_add_executable( NAME BA.05_AMR
        FILES BA.05_AMR.cpp
        DEPENDS hyteg )
if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME BA.05_AMR )
endif()

waLBerla_add_executable( NAME BA.06_ErrorEstimator
        FILES BA.06_ErrorEstimator.cpp
        DEPENDS core hyteg opgen-mass constant_stencil_operator )
if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME BA.06_ErrorEstimator )
endif()
