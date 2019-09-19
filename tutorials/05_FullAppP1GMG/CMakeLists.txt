waLBerla_add_executable( NAME FullAppP1GMG
        FILES P1GMG.cpp
        DEPENDS hyteg core)

waLBerla_link_files_to_builddir( *.prm )
