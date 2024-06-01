add_subdirectory( PlateVelocities ) 
add_subdirectory( SphericalHarmonics )

waLBerla_add_executable( NAME TerraNeo
        FILES TerraNeo.cpp Convection.cpp
        DEPENDS hyteg terraneo core sqlite convection_particles 
        constant_stencil_operator mixed_operator opgen-k_mass opgen-composites-stokes opgen-div_k_grad
        )

waLBerla_link_files_to_builddir( *.prm )
if( HYTEG_BUILD_WITH_ADIOS2 )
    waLBerla_link_files_to_builddir( *.xml )
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()