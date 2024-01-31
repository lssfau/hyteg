waLBerla_add_executable( NAME 08_CahnHilliard
        FILES CahnHilliard.cpp
        DEPENDS hyteg core constantStencilOperator)

waLBerla_link_files_to_builddir( *.prm )
