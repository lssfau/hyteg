if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( MantleConvectionBlendedAnnulus MantleConvectionBlendedAnnulus.cpp )
target_link_libraries( MantleConvectionBlendedAnnulus hyteg convection_particles walberla::sqlite mixed_operator )

add_executable       ( Benchmark_01_CircularAdvection Benchmark_01_CircularAdvection.cpp Helpers.cpp )
target_link_libraries( Benchmark_01_CircularAdvection hyteg convection_particles walberla::sqlite mixed_operator )

add_executable       ( Benchmark_02_SwirlingAdvection Benchmark_02_SwirlingAdvection.cpp Helpers.cpp )
target_link_libraries( Benchmark_02_SwirlingAdvection hyteg convection_particles walberla::sqlite mixed_operator )

add_executable       ( Benchmark_03_BlendedAdvection Benchmark_03_BlendedAdvection.cpp Helpers.cpp )
target_link_libraries( Benchmark_03_BlendedAdvection hyteg convection_particles walberla::sqlite mixed_operator )

add_executable       ( Benchmark_04_BlendedAdvectionDiffusion Benchmark_04_BlendedAdvectionDiffusion.cpp Helpers.cpp )
target_link_libraries( Benchmark_04_BlendedAdvectionDiffusion hyteg convection_particles walberla::sqlite mixed_operator )

add_executable       ( Benchmark_05_PipeScaling Benchmark_05_PipeScaling.cpp Helpers.cpp )
target_link_libraries( Benchmark_05_PipeScaling hyteg convection_particles walberla::sqlite mixed_operator )

if ( HYTEG_BUILD_WITH_PETSC AND HYTEG_TERRANEO_MODULE )
    add_executable       ( MantleConvection MantleConvection.cpp )
    target_link_libraries( MantleConvection hyteg convection_particles walberla::sqlite terraneo mixed_operator )
endif()

add_subdirectory(MantleConvectionRunScripts)
add_subdirectory(MOCScaling)
