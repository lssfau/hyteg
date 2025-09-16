waLBerla_add_executable( NAME MantleConvectionReprod_FrozenVelocity_FGMRES_VCycleBFBT_3D
    FILES MantleConvectionReprod_FrozenVelocity_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MantleConvectionReprod_FrozenVelocity_FGMRES_ScaledMass_3D
    FILES MantleConvectionReprod_FrozenVelocity_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MantleConvectionReprod_FrozenVelocity_FGMRES_WBFBT_3D
    FILES MantleConvectionReprod_FrozenVelocity_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MantleConvectionReprod_FrozenVelocity_FGMRES_diagABFBT_3D
    FILES MantleConvectionReprod_FrozenVelocity_FGMRES_diagABFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)
    
waLBerla_add_executable( NAME MantleConvectionReprod_AsymmetricSystem_FGMRES_VCycleBFBT_3D
    FILES MantleConvectionReprod_AsymmetricSystem_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)
        
waLBerla_add_executable( NAME MantleConvectionReprod_AsymmetricSystem_FGMRES_diagABFBT_3D
    FILES MantleConvectionReprod_AsymmetricSystem_FGMRES_diagABFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)
     
waLBerla_add_executable( NAME MantleConvectionReprod_AsymmetricSystem_FGMRES_ScaledMass_3D
    FILES MantleConvectionReprod_AsymmetricSystem_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)
   
waLBerla_add_executable( NAME MantleConvectionReprod_AsymmetricSystem_FGMRES_WBFBT_3D
    FILES MantleConvectionReprod_AsymmetricSystem_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()