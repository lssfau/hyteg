if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME PlateVelocityDemo
  FILES PlateVelocityDemo.cpp
  DEPENDS hyteg terraneo)

waLBerla_add_executable( NAME PlateAveragingPointsAndWeights
        FILES PlateAveragingPointsAndWeights.cpp
        DEPENDS hyteg terraneo)
