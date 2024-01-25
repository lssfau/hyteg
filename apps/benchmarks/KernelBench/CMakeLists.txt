waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME VertexDoFKernelBench
        FILES VertexDoFKernelBench.cpp
        DEPENDS hyteg core constant_stencil_operator)

waLBerla_add_executable( NAME 3DKernelBench
        FILES 3DKernelBench.cpp
        DEPENDS hyteg core constant_stencil_operator)