if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.xml )

waLBerla_add_executable( NAME GB.01_BlankenbachCase1a
        FILES GB.01_BlankenbachCase1a.cpp
        DEPENDS hyteg terraneo convection_particles sqlite mixed_operator)

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME GB.01_BlankenbachCase1a COMMAND $<TARGET_FILE:GB.01_BlankenbachCase1a> "${CMAKE_BINARY_DIR}/data/param/GB.01_BlankenbachCase1a_CI.prm" PROCESSES 3)
endif()
