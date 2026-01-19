add_executable       ( PrimitiveStorageTransfer PrimitiveStorageTransfer.cpp )
target_link_libraries( PrimitiveStorageTransfer hyteg terraneo walberla::core walberla::sqlite convection_particles opgen-mass )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/config")
endif()

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()
