waLBerla_add_executable( NAME MC_AdvectionDiffusionConvTest
    FILES MC_AdvectionDiffusionConvTest.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MC_PureAdvectionConvTest
    FILES MC_PureAdvectionConvTest.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()