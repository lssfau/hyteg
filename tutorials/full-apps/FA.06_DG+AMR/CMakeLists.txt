if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.ipynb )
waLBerla_link_files_to_builddir( *.prm )

if (HYTEG_BUILD_WITH_PETSC)
    add_executable       ( FA.06_DG+AMR FA.06_DG+AMR.cpp )
    target_link_libraries( FA.06_DG+AMR hyteg walberla::core walberla::sqlite )
    if ( CMAKE_BUILD_TYPE MATCHES "Release" )
        waLBerla_execute_test( NO_MODULE_LABEL NAME FA.06_DG+AMR COMMAND $<TARGET_FILE:FA.06_DG+AMR>
                "${CMAKE_BINARY_DIR}/data/param/FA.06_DG+AMR_CI.prm" )
    endif ()
endif()

