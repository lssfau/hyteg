waLBerla_add_test_executable( P2P1StokesFSGMGSolverTest P2P1StokesFSGMGSolverTest.cpp )
target_link_libraries       ( P2P1StokesFSGMGSolverTest hyteg terraneo walberla::core mixed_operator opgen-k_mass opgen-mass opgen-composites-stokes )
waLBerla_execute_test(NAME P2P1StokesFSGMGSolverTest COMMAND $<TARGET_FILE:P2P1StokesFSGMGSolverTest> LABELS longrun)

waLBerla_add_test_executable( P2P1StokesFSUzawaSolverTest P2P1StokesFSUzawaSolverTest.cpp )
target_link_libraries       ( P2P1StokesFSUzawaSolverTest hyteg terraneo walberla::core mixed_operator opgen-k_mass opgen-mass opgen-composites-stokes )
waLBerla_execute_test(NAME P2P1StokesFSUzawaSolverTest COMMAND $<TARGET_FILE:P2P1StokesFSUzawaSolverTest> LABELS longrun)
