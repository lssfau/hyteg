if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( Poisson-CC-P2-3D Poisson-CC-P2-3D.cpp )
target_link_libraries( Poisson-CC-P2-3D hyteg walberla::core constant_stencil_operator )