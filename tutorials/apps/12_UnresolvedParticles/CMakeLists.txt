if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME 12_UnresolvedParticles
        FILES 12_UnresolvedParticles.cpp
        DEPENDS hyteg core unresolved_particles coupling_hyteg_unresolved_particles )
