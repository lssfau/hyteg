add_executable       ( MantleConvectionReprod_FrozenVelocity_FGMRES_VCycleBFBT_3D MantleConvectionReprod_FrozenVelocity_FGMRES_VCycleBFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_FrozenVelocity_FGMRES_VCycleBFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

add_executable       ( MantleConvectionReprod_FrozenVelocity_FGMRES_ScaledMass_3D MantleConvectionReprod_FrozenVelocity_FGMRES_ScaledMass_3D.cpp )
target_link_libraries( MantleConvectionReprod_FrozenVelocity_FGMRES_ScaledMass_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

add_executable       ( MantleConvectionReprod_FrozenVelocity_FGMRES_WBFBT_3D MantleConvectionReprod_FrozenVelocity_FGMRES_WBFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_FrozenVelocity_FGMRES_WBFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

add_executable       ( MantleConvectionReprod_FrozenVelocity_FGMRES_diagABFBT_3D MantleConvectionReprod_FrozenVelocity_FGMRES_diagABFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_FrozenVelocity_FGMRES_diagABFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )
    
add_executable       ( MantleConvectionReprod_AsymmetricSystem_FGMRES_VCycleBFBT_3D MantleConvectionReprod_AsymmetricSystem_FGMRES_VCycleBFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_AsymmetricSystem_FGMRES_VCycleBFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )
        
add_executable       ( MantleConvectionReprod_AsymmetricSystem_FGMRES_diagABFBT_3D MantleConvectionReprod_AsymmetricSystem_FGMRES_diagABFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_AsymmetricSystem_FGMRES_diagABFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )
     
add_executable       ( MantleConvectionReprod_AsymmetricSystem_FGMRES_ScaledMass_3D MantleConvectionReprod_AsymmetricSystem_FGMRES_ScaledMass_3D.cpp )
target_link_libraries( MantleConvectionReprod_AsymmetricSystem_FGMRES_ScaledMass_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )
   
add_executable       ( MantleConvectionReprod_AsymmetricSystem_FGMRES_WBFBT_3D MantleConvectionReprod_AsymmetricSystem_FGMRES_WBFBT_3D.cpp )
target_link_libraries( MantleConvectionReprod_AsymmetricSystem_FGMRES_WBFBT_3D hyteg walberla::core Boost::headers convection_particles MC-operators elementwise_dof_value_operator )

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()