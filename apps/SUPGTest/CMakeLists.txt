if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *sh )

waLBerla_add_executable( NAME SUPGTestSquare
    FILES  SUPGTestSquare.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles opgen-composites-stokes opgen-div_k_grad opgen-advection opgen-supg_diffusion opgen-supg_advection opgen-full_stokes opgen-terraneo opgen-mass opgen-k_mass opgen-epsilon opgen-divergence opgen-gradient)

waLBerla_add_executable( NAME SUPGTestShell
    FILES  SUPGTestShell.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles opgen-composites-stokes opgen-div_k_grad opgen-supg_mass opgen-advection opgen-supg_diffusion opgen-supg_advection opgen-full_stokes opgen-terraneo opgen-mass opgen-k_mass opgen-epsilon opgen-divergence opgen-gradient)

waLBerla_add_executable( NAME PlumeInRectangle
    FILES  PlumeInRectangle.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles opgen-composites-stokes opgen-div_k_grad opgen-supg_mass opgen-advection opgen-supg_diffusion opgen-supg_advection opgen-full_stokes opgen-terraneo opgen-mass opgen-k_mass opgen-epsilon opgen-divergence opgen-gradient)
