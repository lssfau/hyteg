waLBerla_add_executable( NAME FA.04_IsoviscousConvectionAnnulus
        FILES FA.04_IsoviscousConvectionAnnulus.cpp
        DEPENDS hyteg core sqlite convection_particles mixed_operator)

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME FA.04_IsoviscousConvectionAnnulus )
endif()