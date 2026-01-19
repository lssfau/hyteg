
###################################################################################################
#
# Module to couple unresolved_particles to hyteg
#
###################################################################################################

add_library( coupling_hyteg_unresolved_particles )
target_link_libraries( coupling_hyteg_unresolved_particles PUBLIC walberla::core hyteg unresolved_particles )

add_subdirectory( data )
add_subdirectory( primitivestorage )