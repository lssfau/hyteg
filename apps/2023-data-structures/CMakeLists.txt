if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable(NAME curlCurlConvergence
      FILES curlCurlConvergence.cpp
      DEPENDS hyteg constant_stencil_operator mixed_operator)

waLBerla_add_executable(NAME torusVisualization
      FILES torusVisualization.cpp
      DEPENDS hyteg)
