if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME MantleConvectionBlendedAnnulus
        FILES  MantleConvectionBlendedAnnulus.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

waLBerla_add_executable( NAME Benchmark_01_CircularAdvection
        FILES Benchmark_01_CircularAdvection.cpp Helpers.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

waLBerla_add_executable( NAME Benchmark_02_SwirlingAdvection
        FILES Benchmark_02_SwirlingAdvection.cpp Helpers.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

waLBerla_add_executable( NAME Benchmark_03_BlendedAdvection
        FILES Benchmark_03_BlendedAdvection.cpp Helpers.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

waLBerla_add_executable( NAME Benchmark_04_BlendedAdvectionDiffusion
        FILES Benchmark_04_BlendedAdvectionDiffusion.cpp Helpers.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

waLBerla_add_executable( NAME Benchmark_05_PipeScaling
        FILES Benchmark_05_PipeScaling.cpp Helpers.cpp
        DEPENDS hyteg convection_particles sqlite mixed_operator)

if ( HYTEG_BUILD_WITH_PETSC AND HYTEG_TERRANEO_MODULE )
    waLBerla_add_executable( NAME MantleConvection
            FILES MantleConvection.cpp
            DEPENDS hyteg convection_particles sqlite terraneo mixed_operator)
endif()

add_subdirectory(MantleConvectionRunScripts)
add_subdirectory(MOCScaling)
