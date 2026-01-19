if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( GMRESApp GMRESApp.cpp )
target_link_libraries( GMRESApp hyteg walberla::core constant_stencil_operator )

add_executable       ( StokesGMRES StokesGMRES.cpp )
target_link_libraries( StokesGMRES hyteg mixed_operator )
