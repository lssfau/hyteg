if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
   file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif ()

waLBerla_link_files_to_builddir(*.prm)

waLBerla_add_executable(NAME StokesCubeTransport
      FILES StokesCubeTransport.cpp
      DEPENDS hyteg mixed_operator)
