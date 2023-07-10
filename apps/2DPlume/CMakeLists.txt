if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
   file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME 2DPlume
      FILES 2DPlume.cpp
      DEPENDS hyteg core convection_particles)