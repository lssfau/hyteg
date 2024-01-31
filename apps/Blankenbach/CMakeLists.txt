if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME Blankenbach_Case3
        FILES  Blankenbach_Case3.cpp
        DEPENDS hyteg convection_particles sqlite mixedOperator)
