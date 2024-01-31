if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )


if ( HYTEG_BUILD_WITH_PETSC )
    waLBerla_add_executable( NAME DefectCorrection
            FILES DefectCorrection.cpp
          DEPENDS hyteg core constantStencilOperator)

    waLBerla_add_executable( NAME DefectCorrectionStokes
            FILES DefectCorrectionStokes.cpp
          DEPENDS hyteg core mixedOperator)
endif()
