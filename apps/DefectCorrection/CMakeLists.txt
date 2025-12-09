if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )


if ( HYTEG_BUILD_WITH_PETSC )
    add_executable       ( DefectCorrection DefectCorrection.cpp )
    target_link_libraries( DefectCorrection hyteg walberla::core constant_stencil_operator )

    add_executable       ( DefectCorrectionStokes DefectCorrectionStokes.cpp )
    target_link_libraries( DefectCorrectionStokes hyteg walberla::core mixed_operator )
endif()
