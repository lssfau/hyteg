waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( VertexDoFKernelBench VertexDoFKernelBench.cpp )
target_link_libraries( VertexDoFKernelBench hyteg walberla::core constant_stencil_operator )

add_executable       ( 3DKernelBench 3DKernelBench.cpp )
target_link_libraries( 3DKernelBench hyteg walberla::core constant_stencil_operator )