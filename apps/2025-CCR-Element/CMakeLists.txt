if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME CCRStokes2D
        FILES  CCRStokes2D.cpp
        DEPENDS hyteg core mixed_operator opgen-diffusion manufactured_solutions )
