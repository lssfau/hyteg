if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/benchmark")
endif()

add_compile_options( "-Wno-missing-field-initializers" )

add_subdirectory(operators-used)
waLBerla_link_files_to_builddir(createResidualTablesFromJSON.py)
waLBerla_link_files_to_builddir(count-underflow-from-vtu.py)

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    add_executable       ( benchmarkHockeySticks-residual-convergence benchmarkHockeySticks-residual-convergence.cpp )
    target_link_libraries( benchmarkHockeySticks-residual-convergence walberla::core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    add_executable       ( benchmarkHockeySticks benchmarkHockeySticks.cpp )
    target_link_libraries( benchmarkHockeySticks walberla::core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    add_executable       ( benchmarkHockeySticks-error-threshold benchmarkHockeySticks-error-threshold.cpp )
    target_link_libraries( benchmarkHockeySticks-error-threshold walberla::core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

add_executable       ( benchmarkRuntime benchmarkRuntime.cpp )
target_link_libraries( benchmarkRuntime walberla::core hyteg opgen-mt-diffusion )

add_executable       ( benchmarkRuntimeNoVec benchmarkRuntimeNoVec.cpp )
target_link_libraries( benchmarkRuntimeNoVec walberla::core hyteg opgen-mt-diffusion )

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    add_executable       ( benchmarkFineTuning benchmarkFineTuning.cpp )
    target_link_libraries( benchmarkFineTuning walberla::core hyteg opgen-mt-diffusion )
endif () # WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT

add_executable       ( DoFCounter DoFCounter.cpp )
target_link_libraries( DoFCounter walberla::core hyteg )
