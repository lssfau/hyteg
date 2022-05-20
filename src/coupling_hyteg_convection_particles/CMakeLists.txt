
###################################################################################################
#
# Module to couple convection_particles to hyteg
#
###################################################################################################

add_library( coupling_hyteg_convection_particles )
target_link_libraries( coupling_hyteg_convection_particles PUBLIC core hyteg convection_particles )
target_sources( coupling_hyteg_convection_particles
      PRIVATE
      MMOCTransport.hpp
      )

add_subdirectory( primitivestorage )
add_subdirectory( communication )
