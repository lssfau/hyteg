if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk )
    file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk" )
endif()
if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/csv )
    file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/csv" )
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable(NAME P1GMGSmootherComparison FILES P1GMGSmootherComparison.cpp DEPENDS hyteg core constant_stencil_operator)