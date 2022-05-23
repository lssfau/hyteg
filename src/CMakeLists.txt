add_subdirectory( hyteg )
add_subdirectory( convection_particles )
add_subdirectory( coupling_hyteg_convection_particles )

if( BUILD_TERRANEO_MODULE )
  add_subdirectory( terraneo )
endif()
