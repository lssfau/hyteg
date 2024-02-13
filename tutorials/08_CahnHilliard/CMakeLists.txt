waLBerla_add_executable( NAME 08_CahnHilliard
        FILES CahnHilliard.cpp
        DEPENDS hyteg core constant_stencil_operator)

waLBerla_link_files_to_builddir( *.prm )
