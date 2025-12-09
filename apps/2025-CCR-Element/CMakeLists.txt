if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( CCRStokes2D CCRStokes2D.cpp )
target_link_libraries( CCRStokes2D hyteg walberla::core mixed_operator opgen-diffusion manufactured_solutions )
