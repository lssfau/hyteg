waLBerla_add_executable( NAME PrimitiveStorageTransfer
    FILES PrimitiveStorageTransfer.cpp
    DEPENDS hyteg terraneo core sqlite convection_particles opgen-mass
)

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/config")
endif()

waLBerla_link_files_to_builddir( *.bash )
waLBerla_link_files_to_builddir( *.prm )

if( HYTEG_BUILD_WITH_ADIOS2 )
    waLBerla_link_files_to_builddir( *.xml )
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()
