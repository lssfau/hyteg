if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.pbs )
waLBerla_link_files_to_builddir( *sh )

waLBerla_add_executable( 
    NAME 3DSphShellConvStudies
    FILES  3DSphShellConvStudies.cpp
    DEPENDS hyteg mixed_operator opgen-composites-stokes
)
