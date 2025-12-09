if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( P1CGBenchmark P1CGBenchmark.cpp )
target_link_libraries( P1CGBenchmark hyteg walberla::core constant_stencil_operator )