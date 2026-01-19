add_executable       ( MC_AdvectionDiffusionConvTest MC_AdvectionDiffusionConvTest.cpp )
target_link_libraries( MC_AdvectionDiffusionConvTest hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

add_executable       ( MC_PureAdvectionConvTest MC_PureAdvectionConvTest.cpp )
target_link_libraries( MC_PureAdvectionConvTest hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()