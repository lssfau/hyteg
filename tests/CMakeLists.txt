add_subdirectory( hyteg )
add_subdirectory( convection_particles )
add_subdirectory( manufactured_solutions )
add_subdirectory( unresolved_particles )

if( HYTEG_TERRANEO_MODULE )
  add_subdirectory( terraneo )
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/data/meshes ${CMAKE_CURRENT_BINARY_DIR}/meshes SYMBOLIC)
