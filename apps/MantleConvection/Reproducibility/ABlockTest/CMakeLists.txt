waLBerla_add_executable( NAME MC_ABlockTest_3D
    FILES MC_ABlockTest_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)   

if ( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME MC_ABlockTest_3D_PETSc
    FILES MC_ABlockTest_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)   
endif()

waLBerla_add_executable( NAME MC_ABlockTest_3D_ApproxA
    FILES MC_ABlockTest_3D_ApproxA.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)   

if ( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME MC_ABlockTest_3D_ApproxA_PETSc
    FILES MC_ABlockTest_3D_ApproxA_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)   
endif()

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()