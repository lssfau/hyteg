if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.xml )
waLBerla_link_files_to_builddir( *.prm )

add_executable       ( ADIOS2demo ADIOS2demo.cpp )
target_link_libraries( ADIOS2demo hyteg )

add_executable       ( BenchmarkCheckpointRestore BenchmarkCheckpointRestore.cpp )
target_link_libraries( BenchmarkCheckpointRestore hyteg )
