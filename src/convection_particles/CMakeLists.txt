
###################################################################################################
#
# Module rigid particle dynamics (RPD)
#
###################################################################################################

add_library( convection_particles )
target_link_libraries( convection_particles PUBLIC blockforest core stencil vtk hyteg coupling_hyteg_convection_particles )

add_subdirectory( mpi )
add_subdirectory( common )
add_subdirectory( vtk )
add_subdirectory( sorting )
add_subdirectory( data )
add_subdirectory( domain )
add_subdirectory( collision_detection )
add_subdirectory( kernel )
