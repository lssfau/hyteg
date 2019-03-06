if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )


if ( HHG_BUILD_WITH_PETSC )
    waLBerla_add_executable( NAME DefectCorrection
            FILES DefectCorrection.cpp
            DEPENDS tinyhhg_core core postprocessing)
endif()
