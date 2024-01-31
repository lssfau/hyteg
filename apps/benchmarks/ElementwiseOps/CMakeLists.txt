waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

if( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME PolarLaplacianBenchmark
        FILES PolarLaplacianBenchmark.cpp
        DEPENDS hyteg core)

waLBerla_add_executable( NAME StencilScalingBenchmark
        FILES StencilScalingBenchmark.cpp
        DEPENDS hyteg core constantStencilOperator)
endif()

waLBerla_add_executable( NAME ElementwiseVSConstant
      FILES ElementwiseVSConstant.cpp
      DEPENDS hyteg core constantStencilOperator)
