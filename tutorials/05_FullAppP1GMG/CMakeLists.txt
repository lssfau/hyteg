waLBerla_add_executable( NAME FullAppP1GMG
        FILES P1GMG.cpp
        DEPENDS hyteg core constant_stencil_operator)

waLBerla_link_files_to_builddir( *.prm )
