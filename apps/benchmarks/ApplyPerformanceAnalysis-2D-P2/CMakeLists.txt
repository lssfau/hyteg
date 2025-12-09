if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( ApplyPerformanceAnalysis-2D-P2 ApplyPerformanceAnalysis-2D-P2.cpp AbstractApply.cpp )
target_link_libraries( ApplyPerformanceAnalysis-2D-P2 hyteg walberla::core constant_stencil_operator )