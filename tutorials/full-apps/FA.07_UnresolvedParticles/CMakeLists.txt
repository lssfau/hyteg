if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME FA.07_UnresolvedParticles
        FILES FA.07_UnresolvedParticles.cpp
        DEPENDS hyteg core unresolved_particles coupling_hyteg_unresolved_particles )

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME FA.07_UnresolvedParticles.cpp )
endif()