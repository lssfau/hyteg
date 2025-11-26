if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/checkpoint")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.xml )

waLBerla_add_executable( NAME GB.02_KingCompressible
        FILES GB.02_KingCompressible.cpp
        DEPENDS hyteg terraneo convection_particles sqlite 
        mixed_operator opgen-k_mass opgen-mass opgen-div_k_grad 
        opgen-composites-stokes opgen-shear_heating opgen-terraneo )
