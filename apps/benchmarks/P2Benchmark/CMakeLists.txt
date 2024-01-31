waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME P2Benchmark
        FILES P2Benchmark.cpp
        DEPENDS hyteg core constantStencilOperator)