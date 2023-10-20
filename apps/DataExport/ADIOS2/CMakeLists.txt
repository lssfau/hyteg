if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.xml )
waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME ADIOS2demo
  FILES  ADIOS2demo.cpp
  DEPENDS hyteg )

waLBerla_add_executable( NAME BenchmarkCheckpointRestore
  FILES  BenchmarkCheckpointRestore.cpp
  DEPENDS hyteg )
