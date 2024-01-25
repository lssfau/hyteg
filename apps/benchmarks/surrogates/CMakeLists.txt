waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.job )
waLBerla_link_files_to_builddir( *.sh )

waLBerla_add_executable( NAME SurrogateBenchmarks
        FILES SurrogateBenchmarks.cpp
        DEPENDS hyteg core constant_stencil_operator)

