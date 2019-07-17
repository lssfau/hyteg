waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME P2OperatorBenchmarks
        FILES P2OperatorBenchmarks.cpp
        DEPENDS tinyhhg_core core)
