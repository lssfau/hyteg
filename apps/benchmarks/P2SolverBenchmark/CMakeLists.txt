waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME P2SolverBenchmark
        FILES P2SolverBenchmark.cpp
        DEPENDS hyteg core constant_stencil_operator)