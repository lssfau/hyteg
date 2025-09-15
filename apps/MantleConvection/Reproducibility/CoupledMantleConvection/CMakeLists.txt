waLBerla_add_executable( NAME MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D
    FILES MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D.cpp
    DEPENDS hyteg core Boost::headers sqlite convection_particles constant_stencil_operator mixed_operator MC-composites MC-mass MC-advection MC-diffusion_inv_rho MC-adiabatic_heating MC-composites MC-mass_SUPG MC-advection_supg MC-diffusion_supg MC-adiabatic_heating_supg MC-shear_heating MC-shear_heating_supg MC-full_stokes_scaled_divdiv_vectorial MC-full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base MC-k_mass MC-div_k_grad MC-shear_heating_no_surface MC-shear_heating_supg_no_surface MC-k_p1_coeff_mass elementwise_dof_value_operator MC-operators)   

waLBerla_add_executable( NAME MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_ApproxA
    FILES MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_ApproxA.cpp
    DEPENDS hyteg core Boost::headers sqlite convection_particles constant_stencil_operator mixed_operator MC-composites MC-mass MC-advection MC-diffusion_inv_rho MC-adiabatic_heating MC-composites MC-mass_SUPG MC-advection_supg MC-diffusion_supg MC-adiabatic_heating_supg MC-shear_heating MC-shear_heating_supg MC-full_stokes_scaled_divdiv_vectorial MC-full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base MC-k_mass MC-div_k_grad MC-shear_heating_no_surface MC-shear_heating_supg_no_surface MC-k_p1_coeff_mass MC-operators)   

if ( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_PETSc
    FILES MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers sqlite convection_particles constant_stencil_operator mixed_operator MC-composites MC-mass MC-advection MC-diffusion_inv_rho MC-adiabatic_heating MC-composites MC-mass_SUPG MC-advection_supg MC-diffusion_supg MC-adiabatic_heating_supg MC-shear_heating MC-shear_heating_supg MC-full_stokes_scaled_divdiv_vectorial MC-full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base MC-k_mass MC-div_k_grad MC-shear_heating_no_surface MC-shear_heating_supg_no_surface MC-k_p1_coeff_mass elementwise_dof_value_operator MC-operators)   
endif()

if ( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_ApproxA_PETSc
    FILES MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_3D_ApproxA_PETSc.cpp
    DEPENDS hyteg core Boost::headers sqlite convection_particles constant_stencil_operator mixed_operator MC-composites MC-mass MC-advection MC-diffusion_inv_rho MC-adiabatic_heating MC-composites MC-mass_SUPG MC-advection_supg MC-diffusion_supg MC-adiabatic_heating_supg MC-shear_heating MC-shear_heating_supg MC-full_stokes_scaled_divdiv_vectorial MC-full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base MC-k_mass MC-div_k_grad MC-shear_heating_no_surface MC-shear_heating_supg_no_surface MC-k_p1_coeff_mass MC-operators)   
endif()

if ( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_Arrhenius_3D_PETSc
    FILES MC_TALA_AsymmetricSystem_2ndOrderMMOC_FGMRES_BFBT_HighContrast_Arrhenius_3D_PETSc.cpp
    DEPENDS hyteg core Boost::headers sqlite convection_particles constant_stencil_operator mixed_operator MC-composites MC-mass MC-advection MC-diffusion_inv_rho MC-adiabatic_heating MC-composites MC-mass_SUPG MC-advection_supg MC-diffusion_supg MC-adiabatic_heating_supg MC-shear_heating MC-shear_heating_supg MC-full_stokes_scaled_divdiv_vectorial MC-full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base MC-k_mass MC-div_k_grad MC-shear_heating_no_surface MC-shear_heating_supg_no_surface MC-k_p1_coeff_mass elementwise_dof_value_operator MC-operators)   
endif()

waLBerla_link_files_to_builddir( *.prm )

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()