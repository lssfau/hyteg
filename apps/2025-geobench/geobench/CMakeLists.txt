if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.pbs )
waLBerla_link_files_to_builddir( *.xml )
waLBerla_link_files_to_builddir( *sh )

waLBerla_add_executable( 
    NAME CompressibleTALA
    FILES CompressibleTALA.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles opgen-k_divdiv opgen-mass 
    opgen-k_mass opgen-div_k_grad opgen-epsilon opgen-terraneo opgen-composites-stokes
    opgen-advection
)

waLBerla_add_executable( 
    NAME TosiNonlinear
    FILES  TosiNonlinear.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles 
    opgen-k_mass opgen-composites-stokes
)

waLBerla_add_executable( 
    NAME SphericalShellBench
    FILES  SphericalShellBench.cpp
    DEPENDS hyteg terraneo mixed_operator convection_particles opgen-composites-stokes 
    opgen-div_k_grad opgen-diffusion opgen-advection opgen-supg_mass opgen-supg_diffusion opgen-supg_advection 
    opgen-full_stokes opgen-terraneo opgen-mass opgen-k_mass opgen-epsilon opgen-divergence opgen-gradient
)
