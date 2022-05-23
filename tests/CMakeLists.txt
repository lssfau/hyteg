add_subdirectory( hyteg )
add_subdirectory( convection_particles )

if( BUILD_TERRANEO_MODULE )
  add_subdirectory( terraneo )
endif()
