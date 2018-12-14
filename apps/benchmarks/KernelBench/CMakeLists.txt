waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME VertexDoFKernelBench
        FILES VertexDoFKernelBench.cpp
        DEPENDS tinyhhg_core core)