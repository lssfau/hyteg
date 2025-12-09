if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
   file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( 2DPlume 2DPlume.cpp )
target_link_libraries( 2DPlume hyteg walberla::core mixed_operator convection_particles )