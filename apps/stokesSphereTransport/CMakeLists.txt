if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.sh )

if ( HYTEG_BUILD_WITH_PETSC )

    add_executable       ( StokesSphereTransport StokesSphereTransport.cpp )
    target_link_libraries( StokesSphereTransport hyteg convection_particles mixed_operator )

endif()
