add_executable       ( FA.03_PlumeInCube FA.03_PlumeInCube.cpp )
target_link_libraries( FA.03_PlumeInCube hyteg walberla::core mixed_operator )

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME FA.03_PlumeInCube COMMAND $<TARGET_FILE:FA.03_PlumeInCube>
            "${CMAKE_BINARY_DIR}//data/param/FA.03_PlumeInCube_CI.prm" )
endif()