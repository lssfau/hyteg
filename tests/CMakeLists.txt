add_subdirectory( hyteg )
add_subdirectory( constantStencilOperator )
add_subdirectory( convection_particles )
add_subdirectory( mixedOperator )
add_subdirectory( unresolved_particles )

if( HYTEG_TERRANEO_MODULE )
  add_subdirectory( terraneo )
endif()
