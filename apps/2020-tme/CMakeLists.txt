if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.job )
waLBerla_link_files_to_builddir( *.sh )

add_executable       ( SpMVvsGSBenchmark SpMVvsGSBenchmark.cpp )
target_link_libraries( SpMVvsGSBenchmark hyteg walberla::sqlite mixed_operator )

if (HYTEG_BUILD_WITH_PETSC)
add_executable       ( Benchmark_01_Cube Benchmark_01_Cube.cpp Helpers.cpp )
target_link_libraries( Benchmark_01_Cube hyteg walberla::sqlite mixed_operator )

add_executable       ( Benchmark_02_Y-Pipe Benchmark_02_Y-Pipe.cpp Helpers.cpp )
target_link_libraries( Benchmark_02_Y-Pipe hyteg walberla::sqlite mixed_operator )
endif()
