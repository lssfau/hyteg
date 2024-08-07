if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmark")
endif()

add_compile_options( "-Wno-missing-field-initializers" )

add_subdirectory(operators-used)
waLBerla_link_files_to_builddir(createResidualTablesFromJSON.py)
waLBerla_link_files_to_builddir(count-underflow-from-vtu.py)

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    waLBerla_add_executable( NAME benchmarkHockeySticks-residual-convergence
            FILES benchmarkHockeySticks-residual-convergence.cpp
            DEPENDS core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    waLBerla_add_executable( NAME benchmarkHockeySticks
            FILES benchmarkHockeySticks.cpp
            DEPENDS core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    waLBerla_add_executable( NAME benchmarkHockeySticks-error-threshold
            FILES benchmarkHockeySticks-error-threshold.cpp
            DEPENDS core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

waLBerla_add_executable( NAME benchmarkRuntime
        FILES benchmarkRuntime.cpp
        DEPENDS core hyteg opgen-mt-diffusion )

waLBerla_add_executable( NAME benchmarkRuntimeNoVec
        FILES benchmarkRuntimeNoVec.cpp
        DEPENDS core hyteg opgen-mt-diffusion )

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    waLBerla_add_executable( NAME benchmarkFineTuning
            FILES benchmarkFineTuning.cpp
            DEPENDS core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

waLBerla_add_executable( NAME DoFCounter
        FILES DoFCounter.cpp
        DEPENDS core hyteg )
