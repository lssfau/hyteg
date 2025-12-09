if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( PlateVelocityDemo PlateVelocityDemo.cpp )
target_link_libraries( PlateVelocityDemo hyteg terraneo )

add_executable       ( PlateAveragingPointsAndWeights PlateAveragingPointsAndWeights.cpp )
target_link_libraries( PlateAveragingPointsAndWeights hyteg terraneo )
