waLBerla_add_executable( NAME cg_bubble
        FILES cg_bubble.cpp
        DEPENDS tinyhhg_core core)

waLBerla_add_executable( NAME minres
        FILES minres.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME block_cg
        FILES block_cg.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME stokes_stab
        FILES stokes_stab.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME fullmg
        FILES fullmg.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME cg_P1
        FILES cg_P1.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME cg_P1_varcoeff
        FILES cg_P1_varcoeff.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME cg_P1_tensorvarcoeff
        FILES cg_P1_tensorvarcoeff.cpp
        DEPENDS tinyhhg_core)

waLBerla_execute_test(NO_MODULE_LABEL NAME cg_bubble       CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME minres          CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME block_cg        CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME stokes_stab     CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME fullmg          CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME cg_P1           CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME cg_P1_varcoeff  CONFIGURATIONS Release RelWithDbgInfo)

waLBerla_add_executable( NAME chorin
        FILES chorin_test.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME stokes_mini
        FILES stokes_mini.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME stokes_stab_varcoeff
        FILES stokes_stab_varcoeff.cpp
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

waLBerla_add_executable( NAME cg_P1_varcoords
        FILES cg_P1_varcoords.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME cg_P2
        FILES cg_P2.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME gmg_P2
        FILES gmg_P2.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME GeometryBlending
        FILES GeometryBlending.cpp
        DEPENDS tinyhhg_core)

if( HHG_BUILD_WITH_EIGEN )
waLBerla_add_executable( NAME cg_P1_varcoeff_poly
        FILES cg_P1_varcoeff_poly.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME gmg_varcoeff_vs_poly
        FILES gmg_varcoeff_vs_poly.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME gmg_blending
        FILES gmg_blending.cpp
        DEPENDS tinyhhg_core)
endif()


if( HHG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME stokes_mini_petsc
        FILES stokes_mini_petsc.cpp
        DEPENDS tinyhhg_core)
endif()

waLBerla_add_executable( NAME CompareInterpolate
        FILES benchmarks/CompareInterpolate.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME ApplyBenchmark
        FILES benchmarks/ApplyBenchmark.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME IndexFunctionVSManualVersion
        FILES benchmarks/IndexFunctionVSManualVersion.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME show_mesh
        FILES show_mesh.cpp
        DEPENDS tinyhhg_core)

waLBerla_add_executable( NAME StokesP2P1
        FILES StokesP2P1.cpp
        DEPENDS tinyhhg_core)


#add_executable (tinytest_fmg fmg_test.cpp)
#target_link_libraries (tinytest_fmg tinyhhg ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})
#
#add_executable (tinytest_chorin chorin_test.cpp)
#target_link_libraries (tinytest_chorin tinyhhg ${MPI_CXX_LIBRARIES})
