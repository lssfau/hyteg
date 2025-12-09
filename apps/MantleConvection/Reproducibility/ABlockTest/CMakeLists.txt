add_executable       ( MC_ABlockTest_3D MC_ABlockTest_3D.cpp )
target_link_libraries( MC_ABlockTest_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )   

if ( HYTEG_BUILD_WITH_PETSC )
add_executable       ( MC_ABlockTest_3D_PETSc MC_ABlockTest_3D_PETSc.cpp )
target_link_libraries( MC_ABlockTest_3D_PETSc hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )   
endif()

add_executable       ( MC_ABlockTest_3D_ApproxA MC_ABlockTest_3D_ApproxA.cpp )
target_link_libraries( MC_ABlockTest_3D_ApproxA hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )   

if ( HYTEG_BUILD_WITH_PETSC )
add_executable       ( MC_ABlockTest_3D_ApproxA_PETSc MC_ABlockTest_3D_ApproxA_PETSc.cpp )
target_link_libraries( MC_ABlockTest_3D_ApproxA_PETSc hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )   
endif()

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()