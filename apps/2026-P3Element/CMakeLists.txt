if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( P3ShapeFunctionDemo P3ShapeFunctionDemo.cpp )
target_link_libraries( P3ShapeFunctionDemo hyteg walberla::core )
