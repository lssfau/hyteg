waLBerla_add_executable( NAME IsoviscousConvectionAnnulus
        FILES IsoviscousConvectionAnnulus.cpp
        DEPENDS hyteg core sqlite convection_particles mixed_operator)

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
