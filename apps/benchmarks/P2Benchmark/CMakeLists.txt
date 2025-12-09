waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P2Benchmark P2Benchmark.cpp )
target_link_libraries( P2Benchmark hyteg walberla::core constant_stencil_operator )