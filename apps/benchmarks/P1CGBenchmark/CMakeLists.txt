waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME P1CGBenchmark
        FILES P1CGBenchmark.cpp
        DEPENDS tinyhhg_core core)