add_compile_options( "-std=c++20" )
add_flag( CMAKE_CXX_FLAGS "-Wno-ambiguous-reversed-operator" )

add_executable       ( TerraNeo Convection.cpp )
target_link_libraries( TerraNeo
    hyteg terraneo walberla::core walberla::sqlite convection_particles
    constant_stencil_operator mixed_operator opgen-mass opgen-k_mass
    opgen-composites-stokes opgen-div_k_grad opgen-terraneo
)

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/config")
endif()
# waLBerla_link_files_to_builddir( config/*.prm )
waLBerla_link_files_to_builddir( *.bash )
waLBerla_link_files_to_builddir( *.prm )
if( HYTEG_BUILD_WITH_ADIOS2 )
    waLBerla_link_files_to_builddir( *.xml )
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()
