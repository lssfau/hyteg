if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.job )
waLBerla_link_files_to_builddir( *.sh )

waLBerla_add_executable( NAME flexcycles
        FILES flexcycles.cpp
        DEPENDS hyteg core opgen-div_k_grad constant_stencil_operator)
 