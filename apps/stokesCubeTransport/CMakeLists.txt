if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
   file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif ()

waLBerla_link_files_to_builddir(*.prm)

add_executable       ( StokesCubeTransport StokesCubeTransport.cpp )
target_link_libraries( StokesCubeTransport hyteg mixed_operator )
