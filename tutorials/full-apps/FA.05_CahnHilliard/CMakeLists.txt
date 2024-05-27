waLBerla_add_executable( NAME FA.05_CahnHilliard
        FILES FA.05_CahnHilliard.cpp
        DEPENDS hyteg core constant_stencil_operator)

waLBerla_link_files_to_builddir( *.prm )
