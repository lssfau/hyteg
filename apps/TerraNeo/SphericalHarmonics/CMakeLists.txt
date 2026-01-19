if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( SPHdemo SPHdemo.cpp )
target_link_libraries( SPHdemo hyteg terraneo )

if (NOT APPLE)
  add_executable       ( EigFunLaplace EigFunLaplace.cpp )
  target_link_libraries( EigFunLaplace hyteg terraneo )
endif ()

add_executable       ( bmark bmark.cpp )
target_link_libraries( bmark hyteg terraneo mixed_operator )
