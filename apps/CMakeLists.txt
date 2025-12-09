add_executable       ( stokes_flow_through_porous_med stokes_flow_through_porous_med.cpp )
target_link_libraries( stokes_flow_through_porous_med hyteg mixed_operator )

add_executable       ( StokesFlowSolverComparison StokesFlowSolverComparison.cpp )
target_link_libraries( StokesFlowSolverComparison hyteg mixed_operator )

add_executable       ( fullmg fullmg.cpp )
target_link_libraries( fullmg hyteg mixed_operator )

add_executable       ( Chorin Chorin.cpp )
target_link_libraries( Chorin hyteg mixed_operator )

add_executable       ( p0p1Transport p0p1Transport.cpp )
target_link_libraries( p0p1Transport hyteg )

add_executable       ( stokes_stab_uzawa stokes_stab_uzawa.cpp )
target_link_libraries( stokes_stab_uzawa hyteg mixed_operator )

add_executable       ( geophysics_plume geophysics_plume.cpp )
target_link_libraries( geophysics_plume hyteg mixed_operator )

add_executable       ( GeometryBlending GeometryBlending.cpp )
target_link_libraries( GeometryBlending hyteg mixed_operator )

add_executable       ( show_mesh show_mesh.cpp )
target_link_libraries( show_mesh hyteg )

add_executable       ( PolarLaplacian PolarLaplacian.cpp )
target_link_libraries( PolarLaplacian hyteg )

add_executable       ( PrimitiveSizes PrimitiveSizes.cpp )
target_link_libraries( PrimitiveSizes hyteg )

add_executable       ( refCellCouplingCount refCellCouplingCount.cpp )
target_link_libraries( refCellCouplingCount hyteg constant_stencil_operator )

add_executable       ( cg_P1_blending cg_P1_blending.cpp )
target_link_libraries( cg_P1_blending hyteg )

add_executable       ( gmg_blending gmg_blending.cpp )
target_link_libraries( gmg_blending hyteg )

if( HYTEG_BUILD_WITH_PETSC )
add_executable       ( exportOperatorMatrix exportOperatorMatrix.cpp )
target_link_libraries( exportOperatorMatrix hyteg mixed_operator opgen-diffusion )

add_executable       ( annulus annulus.cpp )
target_link_libraries( annulus hyteg )
endif()

add_subdirectory(2020-moc)
add_subdirectory(2020-scaling-workshop)
add_subdirectory(2020-tme)
add_subdirectory(2022-eg-varvisc)
add_subdirectory(2023-bauer-mt)
add_subdirectory(2023-data-structures)
add_subdirectory(2023-zikeli-mt)
add_subdirectory(2024-parthasarathy-mgoptimisation)
add_subdirectory(2DPlume)
add_subdirectory(Blankenbach)
add_subdirectory(ConvectionCell2D)
add_subdirectory(DataExport)
add_subdirectory(DefectCorrection)
add_subdirectory(GKBTestProblemGeneration)
add_subdirectory(GMRESTestApp)
add_subdirectory(MG_with_mixed_BCs)
add_subdirectory(MultigridStudies)
add_subdirectory(Surrogates)
add_subdirectory(AMR)
add_subdirectory(benchmarks)
add_subdirectory(blockOperators)
add_subdirectory(profiling)
add_subdirectory(scaling)
add_subdirectory(smootherPerformance)
add_subdirectory(stokesCubeTransport)
add_subdirectory(stokesFreeSlip)
add_subdirectory(stokesSphere)
add_subdirectory(stokesSphereTransport)
add_subdirectory(surfaces)

if( HYTEG_BUILD_WITH_PETSC )
  add_subdirectory(2021-tokamak)
  add_subdirectory(2025-CCR-Element)
endif()

if( HYTEG_TERRANEO_MODULE )
  add_subdirectory(TerraNeo)
endif()

if( HYTEG_MANTLECONVECTION_APP )
  add_subdirectory( MantleConvection )
endif()

if( HYTEG_BUILD_WITH_PYTHON3 )
  add_subdirectory(2024-convbench)
endif()
