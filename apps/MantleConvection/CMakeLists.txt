add_subdirectory(Utility)
add_subdirectory(Reproducibility)

# Main models

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_Arrhenius_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_Arrhenius_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_Arrhenius_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_Arrhenius_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_Arrhenius_3D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_Arrhenius_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

if ( HYTEG_BUILD_WITH_PETSC )

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_Arrhenius_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_Arrhenius_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_Arrhenius_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_Arrhenius_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_Arrhenius_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_Arrhenius_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

endif()


# Rest of the models that serve as building blocks and are not compiled during CI to avoid
# large build times. Since HYTEG_DONT_BUILD_DURING_CI is ON by default, the build targets
# exists during normal usage unless you manually change HYTEG_DONT_BUILD_DURING_CI to OFF.

if ( HYTEG_DONT_BUILD_DURING_CI )

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_3D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_3D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_3D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_3D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_3D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_3D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_3D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_3D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_3D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_3D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_2D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_2D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_2D
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_2D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_2D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_2D
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_2D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_2D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_2D
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_2D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_2D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_2D
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_2D.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

if ( HYTEG_BUILD_WITH_PETSC )

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_3D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_3D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_ScaledMass_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_VCycleBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_MMOC_FGMRES_WBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_ScaledMass_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_VCycleBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_2D_PETSc
    FILES MCMAIN_TALA_AsymmetricSystem_SUPG_FGMRES_WBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_ScaledMass_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_VCycleBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_MMOC_FGMRES_WBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_ScaledMass_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_VCycleBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

waLBerla_add_executable( NAME MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_2D_PETSc
    FILES MCMAIN_TALA_FrozenVelocity_SUPG_FGMRES_WBFBT_2D_PETSc.cpp
    DEPENDS hyteg core Boost::headers convection_particles MC-operators elementwise_dof_value_operator)

endif()

endif()

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()