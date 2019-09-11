if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME TauExtrapolation3D
        FILES TauExtrapolation3D.cpp
        DEPENDS hyteg core sqlite)

