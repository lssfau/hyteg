waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.sh )

waLBerla_add_executable( NAME SnoopFilterIssueBenchmark
        FILES SnoopFilterIssueBenchmark.cpp
        DEPENDS hyteg core sqlite constantStencilOperator)