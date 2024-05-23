waLBerla_add_executable( NAME FA.03_PlumeInCube
        FILES FA.03_PlumeInCube.cpp
        DEPENDS hyteg core mixed_operator)

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
