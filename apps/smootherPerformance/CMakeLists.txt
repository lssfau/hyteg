if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk )
    file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk" )
endif()
if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/csv )
    file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/csv" )
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P1GMGSmootherComparison P1GMGSmootherComparison.cpp )
target_link_libraries( P1GMGSmootherComparison hyteg walberla::core constant_stencil_operator )