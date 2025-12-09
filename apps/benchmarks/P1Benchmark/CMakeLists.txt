waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P1Benchmark P1Benchmark.cpp )
target_link_libraries( P1Benchmark hyteg walberla::core constant_stencil_operator )