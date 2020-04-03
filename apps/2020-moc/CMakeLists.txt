if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME MantleConvectionBlendedAnnulus
        FILES  MantleConvectionBlendedAnnulus.cpp
        DEPENDS hyteg convection_particles)
