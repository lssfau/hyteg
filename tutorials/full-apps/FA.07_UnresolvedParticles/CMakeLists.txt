if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( FA.07_UnresolvedParticles FA.07_UnresolvedParticles.cpp )
target_link_libraries( FA.07_UnresolvedParticles hyteg walberla::core unresolved_particles coupling_hyteg_unresolved_particles )

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME FA.07_UnresolvedParticles.cpp COMMAND $<TARGET_FILE:FA.07_UnresolvedParticles>
            "${CMAKE_BINARY_DIR}/data/param/FA.07_UnresolvedParticles_CI.prm" )
endif()