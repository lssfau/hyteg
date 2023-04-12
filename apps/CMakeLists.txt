waLBerla_add_executable( NAME stokes_flow_through_porous_med
        FILES stokes_flow_through_porous_med.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME StokesFlowSolverComparison
        FILES StokesFlowSolverComparison.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME fullmg
        FILES fullmg.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME Chorin
        FILES Chorin.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME p0p1Transport
        FILES p0p1Transport.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME stokes_stab_uzawa
        FILES stokes_stab_uzawa.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME geophysics_plume
        FILES geophysics_plume.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME GeometryBlending
        FILES GeometryBlending.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME show_mesh
        FILES show_mesh.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME PolarLaplacian
        FILES  PolarLaplacian.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME PrimitiveSizes
        FILES PrimitiveSizes.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME refCellCouplingCount
        FILES refCellCouplingCount.cpp
        DEPENDS hyteg)

waLBerla_add_executable(NAME cg_P1_blending
      FILES cg_P1_blending.cpp
      DEPENDS hyteg)

waLBerla_add_executable(NAME gmg_blending
      FILES gmg_blending.cpp
      DEPENDS hyteg)

if (NOT APPLE)
   # apple clang does not contain expint which is used in this app
   waLBerla_add_executable(NAME adaptiveRefinement
         FILES adaptiveRefinement.cpp
         DEPENDS hyteg)
endif ()

if( HYTEG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME exportOperatorMatrix
        FILES exportOperatorMatrix.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME annulus
        FILES annulus.cpp
        DEPENDS hyteg)
endif()

add_subdirectory(benchmarks)
add_subdirectory(Blankenbach)
add_subdirectory(MG_with_mixed_BCs)
add_subdirectory(scaling)
add_subdirectory(stokesSphere)
add_subdirectory(stokesCubeTransport)
add_subdirectory(stokesFreeSlip)
add_subdirectory(stokesSphereTransport)
add_subdirectory(MultigridStudies)
add_subdirectory(DefectCorrection)
add_subdirectory(smootherPerformance)
add_subdirectory(ConvectionCell2D)
add_subdirectory(2020-moc)
add_subdirectory(2020-tme)
add_subdirectory(2020-scaling-workshop)
add_subdirectory(GMRESTestApp)
add_subdirectory(GKBTestProblemGeneration)

if( HYTEG_BUILD_WITH_PETSC )
  add_subdirectory(2021-tokamak)
endif()

if( HYTEG_TERRANEO_MODULE )
  add_subdirectory(TerraNeo)
endif()

add_subdirectory(surfaces)

add_subdirectory(Surrogates)

add_subdirectory(blockOperators)

add_subdirectory(2022-eg-varvisc)
