waLBerla_add_executable( NAME stokes_flow_through_porous_med
        FILES stokes_flow_through_porous_med.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME StokesFlowSolverComparison
        FILES StokesFlowSolverComparison.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME fullmg
        FILES fullmg.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME Chorin
        FILES Chorin.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME dg0_transport
        FILES dg0_transport.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME stokes_stab_uzawa
        FILES stokes_stab_uzawa.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME geophysics_plume
        FILES geophysics_plume.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME GeometryBlending
        FILES GeometryBlending.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME show_mesh
        FILES show_mesh.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME PolarLaplacian
        FILES  polar_laplacian.cpp
        DEPENDS tinyhhg_core)

if( HHG_BUILD_WITH_EIGEN )
    waLBerla_add_executable( NAME cg_P1_blending
            FILES cg_P1_blending.cpp
            DEPENDS tinyhhg_core)

    waLBerla_add_executable( NAME gmg_blending
            FILES gmg_blending.cpp
            DEPENDS tinyhhg_core)
endif()

add_subdirectory(benchmarks)
add_subdirectory(stokesSphere)
add_subdirectory(stokesCubeTransport)
