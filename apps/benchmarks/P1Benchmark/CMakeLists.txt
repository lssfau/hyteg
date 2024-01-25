waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME P1Benchmark
        FILES P1Benchmark.cpp
        DEPENDS hyteg core constant_stencil_operator)