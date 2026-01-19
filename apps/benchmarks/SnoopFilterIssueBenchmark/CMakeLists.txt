waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.sh )

add_executable       ( SnoopFilterIssueBenchmark SnoopFilterIssueBenchmark.cpp )
target_link_libraries( SnoopFilterIssueBenchmark hyteg walberla::core walberla::sqlite constant_stencil_operator )