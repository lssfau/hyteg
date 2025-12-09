waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

if( HYTEG_BUILD_WITH_PETSC )
add_executable       ( PolarLaplacianBenchmark PolarLaplacianBenchmark.cpp )
target_link_libraries( PolarLaplacianBenchmark hyteg walberla::core )

add_executable       ( StencilScalingBenchmark StencilScalingBenchmark.cpp )
target_link_libraries( StencilScalingBenchmark hyteg walberla::core constant_stencil_operator )
endif()

add_executable       ( ElementwiseVSConstant ElementwiseVSConstant.cpp )
target_link_libraries( ElementwiseVSConstant hyteg walberla::core constant_stencil_operator )
