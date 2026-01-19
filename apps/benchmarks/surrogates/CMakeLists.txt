waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.job )
waLBerla_link_files_to_builddir( *.sh )

add_executable       ( SurrogateBenchmarks SurrogateBenchmarks.cpp )
target_link_libraries( SurrogateBenchmarks hyteg walberla::core constant_stencil_operator )

