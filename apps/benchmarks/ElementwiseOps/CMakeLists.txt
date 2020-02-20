waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

if( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME PolarLaplacianBenchmark
        FILES PolarLaplacianBenchmark.cpp
        DEPENDS hyteg core)
endif()
