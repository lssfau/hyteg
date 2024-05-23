waLBerla_add_executable( NAME FA.01_GeometricMultigrid
        FILES FA.01_GeometricMultigrid.cpp
        DEPENDS hyteg core constant_stencil_operator)

waLBerla_link_files_to_builddir( *.prm )
