waLBerla_add_executable( NAME cg_bubble
        FILES cg_bubble.cpp
        DEPENDS tinyhhg_core fmt core)

waLBerla_add_executable( NAME minres
        FILES minres.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME block_cg
        FILES block_cg.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME stokes_stab
        FILES stokes_stab.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME fullmg
        FILES fullmg.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME cg_P1
        FILES cg_P1.cpp
        DEPENDS tinyhhg_core fmt)



waLBerla_execute_test(NO_MODULE_LABEL NAME cg_bubble   CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME minres      CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME block_cg    CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME stokes_stab CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME fullmg      CONFIGURATIONS Release RelWithDbgInfo)
waLBerla_execute_test(NO_MODULE_LABEL NAME cg_P1       CONFIGURATIONS Release RelWithDbgInfo)

#waLBerla_add_executable( NAME tinytest_chorin
#        FILES chorin_test.cpp
#        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME stokes_mini
        FILES stokes_mini.cpp
        DEPENDS tinyhhg_core fmt)

if( HHG_BUILD_WITH_PETSC )
waLBerla_add_executable( NAME exact_solver
        FILES exact_solver.cpp
        DEPENDS tinyhhg_core fmt)
endif()

#add_executable (tinytest_fmg fmg_test.cpp)
#target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})
#
#add_executable (tinytest_chorin chorin_test.cpp)
#target_link_libraries (tinytest_chorin tinyhhg fmt ${MPI_CXX_LIBRARIES})
