if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( ConvectionCell2D ConvectionCell2D.cpp )
target_link_libraries( ConvectionCell2D hyteg convection_particles mixed_operator )
