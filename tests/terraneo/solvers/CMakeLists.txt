waLBerla_compile_test(FILES P2P1StokesFSGMGSolverTest.cpp DEPENDS hyteg terraneo core mixed_operator opgen-k_mass opgen-mass opgen-composites-stokes)
waLBerla_execute_test(NAME P2P1StokesFSGMGSolverTest COMMAND $<TARGET_FILE:P2P1StokesFSGMGSolverTest> LABELS longrun)

waLBerla_compile_test(FILES P2P1StokesFSUzawaSolverTest.cpp DEPENDS hyteg terraneo core mixed_operator opgen-k_mass opgen-mass opgen-composites-stokes)
waLBerla_execute_test(NAME P2P1StokesFSUzawaSolverTest COMMAND $<TARGET_FILE:P2P1StokesFSUzawaSolverTest> LABELS longrun)
