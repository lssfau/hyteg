if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/checkpoint)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/checkpoint")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.xml )

add_executable       ( GB.02_KingCompressible GB.02_KingCompressible.cpp )
target_link_libraries( GB.02_KingCompressible hyteg terraneo convection_particles walberla::sqlite mixed_operator
                                              opgen-k_mass opgen-mass opgen-div_k_grad opgen-composites-stokes 
                                              opgen-shear_heating opgen-terraneo )

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME GB.02_KingCompressible COMMAND $<TARGET_FILE:GB.02_KingCompressible> "${CMAKE_BINARY_DIR}/data/param/GB.02_KingCompressible_CI.prm" PROCESSES 4)
endif()
