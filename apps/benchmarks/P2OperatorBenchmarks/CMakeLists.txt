waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P2OperatorBenchmarks P2OperatorBenchmarks.cpp )
target_link_libraries( P2OperatorBenchmarks hyteg walberla::core constant_stencil_operator )

add_subdirectory(2019_tme)