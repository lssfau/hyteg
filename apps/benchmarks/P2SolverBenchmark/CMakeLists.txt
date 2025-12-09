waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P2SolverBenchmark P2SolverBenchmark.cpp )
target_link_libraries( P2SolverBenchmark hyteg walberla::core constant_stencil_operator )