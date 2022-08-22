if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.ipynb )

waLBerla_add_executable( NAME 10_DGAMR
        FILES 10_DGAMR.cpp
        DEPENDS hyteg core sqlite)
