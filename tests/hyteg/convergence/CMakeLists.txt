waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_test_executable( ElementwiseDivKGradCGConvergenceTest ElementwiseDivKGradCGConvergenceTest.cpp )
target_link_libraries       ( ElementwiseDivKGradCGConvergenceTest hyteg walberla::core opgen-div_k_grad )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradCGConvergenceTest)
endif ()

waLBerla_add_test_executable( ElementwiseDivKGradNeumannCGConvergenceTest ElementwiseDivKGradNeumannCGConvergenceTest.cpp )
target_link_libraries       ( ElementwiseDivKGradNeumannCGConvergenceTest hyteg walberla::core opgen-div_k_grad )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradNeumannCGConvergenceTest)
endif ()

waLBerla_add_test_executable( P2P1StokesPreconditionersConvergenceTest P2P1StokesPreconditionersConvergenceTest.cpp )
target_link_libraries       ( P2P1StokesPreconditionersConvergenceTest hyteg walberla::core opgen-k_mass opgen-composites-stokes constant_stencil_operator )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME P2P1StokesPreconditionersConvergenceTest)
    waLBerla_execute_test(NAME P2P1StokesPreconditionersConvergenceTestLongrun COMMAND $<TARGET_FILE:P2P1StokesPreconditionersConvergenceTest> --longrun LABELS longrun)
endif ()

waLBerla_add_test_executable( P1VectorLaplaceCGTest P1VectorLaplaceCGTest.cpp )
target_link_libraries       ( P1VectorLaplaceCGTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P1VectorLaplaceCGTest)

waLBerla_add_test_executable( P1StokesMinResConvergenceTest P1StokesMinResConvergenceTest.cpp )
target_link_libraries       ( P1StokesMinResConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P1StokesMinResConvergenceTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P1P1StokesUzawaConvergenceTest P1P1StokesUzawaConvergenceTest.cpp )
    target_link_libraries       ( P1P1StokesUzawaConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P1P1StokesUzawaConvergenceTest)
endif ()

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_add_test_executable( P1Stokes3DMinResConvergenceTest P1Stokes3DMinResConvergenceTest.cpp )
    target_link_libraries       ( P1Stokes3DMinResConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P1Stokes3DMinResConvergenceTest)

    waLBerla_add_test_executable( P2P1Stokes3DMinResConvergenceTest P2P1Stokes3DMinResConvergenceTest.cpp )
    target_link_libraries       ( P2P1Stokes3DMinResConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P2P1Stokes3DMinResConvergenceTest)

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_add_test_executable( P1P1Stokes3DUzawaConvergenceTest P1P1Stokes3DUzawaConvergenceTest.cpp )
        target_link_libraries       ( P1P1Stokes3DUzawaConvergenceTest hyteg walberla::core mixed_operator )
        waLBerla_execute_test(NAME P1P1Stokes3DUzawaConvergenceTest)
    endif ()

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_add_test_executable( P2P1Stokes3DUzawaConvergenceTest P2P1Stokes3DUzawaConvergenceTest.cpp )
        target_link_libraries       ( P2P1Stokes3DUzawaConvergenceTest hyteg walberla::core mixed_operator )
        waLBerla_execute_test(NAME P2P1Stokes3DUzawaConvergenceTest)
    endif ()
endif ()

waLBerla_add_test_executable( P2ElementwiseCGConvergenceTest P2ElementwiseCGConvergenceTest.cpp )
target_link_libraries       ( P2ElementwiseCGConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P2ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( ElementwiseEpsilonMinResConvergenceTest ElementwiseEpsilonMinResConvergenceTest.cpp )
target_link_libraries       ( ElementwiseEpsilonMinResConvergenceTest hyteg walberla::core mixed_operator opgen-composites-stokes )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTest)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTestLongrun COMMAND $<TARGET_FILE:ElementwiseEpsilonMinResConvergenceTest> --longrun LABELS longrun)
    endif ()
endif ()

waLBerla_add_test_executable( P2P1StokesMinResConvergenceTest P2P1StokesMinResConvergenceTest.cpp )
target_link_libraries       ( P2P1StokesMinResConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P2P1StokesMinResConvergenceTest)

waLBerla_add_test_executable( P2P1StokesSchurCGConvergenceTest P2P1StokesSchurCGConvergenceTest.cpp )
target_link_libraries       ( P2P1StokesSchurCGConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P2P1StokesSchurCGConvergenceTest)

waLBerla_add_test_executable( P2P1ElementwiseUzawaConvergenceTest P2P1ElementwiseUzawaConvergenceTest.cpp )
target_link_libraries       ( P2P1ElementwiseUzawaConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P2P1ElementwiseUzawaConvergenceTest)

waLBerla_add_test_executable( P2P1UzawaConvergenceTest P2P1UzawaConvergenceTest.cpp )
target_link_libraries       ( P2P1UzawaConvergenceTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME P2P1UzawaConvergenceTest)

waLBerla_add_test_executable( StokesVarViscPopovTest StokesVarViscPopovTest.cpp )
target_link_libraries       ( StokesVarViscPopovTest hyteg walberla::core mixed_operator )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME StokesVarViscPopovTest)
    endif ()
endif ()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_add_test_executable( P2P1Stokes2DGKBConvergenceTest P2P1Stokes2DGKBConvergenceTest.cpp )
    target_link_libraries       ( P2P1Stokes2DGKBConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P2P1Stokes2DGKBConvergenceTest)
endif()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_add_test_executable( P2P1ChannelGKBConvergenceTest P2P1ChannelGKBConvergenceTest.cpp )
    target_link_libraries       ( P2P1ChannelGKBConvergenceTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P2P1ChannelGKBConvergenceTest)
endif()

waLBerla_add_test_executable( P1CGConvergenceTest P1CGConvergenceTest.cpp )
target_link_libraries       ( P1CGConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1CGConvergenceTest)
waLBerla_execute_test(NAME P1CGConvergenceTestMPI COMMAND $<TARGET_FILE:P1CGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1CG3DConvergenceTest P1CG3DConvergenceTest.cpp )
target_link_libraries       ( P1CG3DConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1CG3DConvergenceTest)
waLBerla_execute_test(NAME P1CG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1CG3DConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1GMGConvergenceTest P1GMGConvergenceTest.cpp )
target_link_libraries       ( P1GMGConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1GMGConvergenceTest)
waLBerla_execute_test(NAME P1GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1FASConvergenceTest P1FASConvergenceTest.cpp )
target_link_libraries       ( P1FASConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1FASConvergenceTest)
waLBerla_execute_test(NAME P1FASConvergenceTestMPI COMMAND $<TARGET_FILE:P1FASConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1TauExtrapolationConvergenceTest P1TauExtrapolationConvergenceTest.cpp )
target_link_libraries       ( P1TauExtrapolationConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTest)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTestMPI COMMAND $<TARGET_FILE:P1TauExtrapolationConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1GMG3DConvergenceTest P1GMG3DConvergenceTest.cpp )
target_link_libraries       ( P1GMG3DConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1GMG3DConvergenceTest)
waLBerla_execute_test(NAME P1GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMG3DConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P2GMGConvergenceTest P2GMGConvergenceTest.cpp )
target_link_libraries       ( P2GMGConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2GMGConvergenceTest)
waLBerla_execute_test(NAME P2GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P2GMG3DConvergenceTest P2GMG3DConvergenceTest.cpp )
target_link_libraries       ( P2GMG3DConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2GMG3DConvergenceTest)
waLBerla_execute_test(NAME P2GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMG3DConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_add_test_executable( P1GSConvergenceTest P1GSConvergenceTest.cpp )
    target_link_libraries       ( P1GSConvergenceTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P1GSConvergenceTest)
endif()

waLBerla_add_test_executable( P1JacobiConvergenceTest P1JacobiConvergenceTest.cpp )
target_link_libraries       ( P1JacobiConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1JacobiConvergenceTest)

waLBerla_add_test_executable( P1MinResConvergenceTest P1MinResConvergenceTest.cpp )
target_link_libraries       ( P1MinResConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1MinResConvergenceTest)

waLBerla_add_test_executable( P2CGConvergenceTest P2CGConvergenceTest.cpp )
target_link_libraries       ( P2CGConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2CGConvergenceTest)
waLBerla_execute_test(NAME P2CGConvergenceTestMPI COMMAND $<TARGET_FILE:P2CGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P1ElementwiseCGConvergenceTest P1ElementwiseCGConvergenceTest.cpp )
target_link_libraries       ( P1ElementwiseCGConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P1ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P2P1FGMRESStokesTest P2P1FGMRESStokesTest.cpp )
target_link_libraries       ( P2P1FGMRESStokesTest hyteg walberla::core constant_stencil_operator mixed_operator )
waLBerla_execute_test(NAME P2P1FGMRESStokesTest)

waLBerla_add_test_executable( P2P1GMRESStokesTest P2P1GMRESStokesTest.cpp )
target_link_libraries       ( P2P1GMRESStokesTest hyteg walberla::core constant_stencil_operator mixed_operator )
waLBerla_execute_test(NAME P2P1GMRESStokesTest)

waLBerla_add_test_executable( P2GMRESConvergenceTest P2GMRESConvergenceTest.cpp )
target_link_libraries       ( P2GMRESConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2GMRESConvergenceTest)

waLBerla_add_test_executable( P2FGMRESConvergenceTest P2FGMRESConvergenceTest.cpp )
target_link_libraries       ( P2FGMRESConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2FGMRESConvergenceTest)

waLBerla_add_test_executable( P1GMRESConvergenceTest P1GMRESConvergenceTest.cpp )
target_link_libraries       ( P1GMRESConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1GMRESConvergenceTest)

waLBerla_add_test_executable( P1FGMRESConvergenceTest P1FGMRESConvergenceTest.cpp )
target_link_libraries       ( P1FGMRESConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P1FGMRESConvergenceTest)

waLBerla_add_test_executable( P2SmoothConvergenceTest P2SmoothConvergenceTest.cpp )
target_link_libraries       ( P2SmoothConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2SmoothConvergenceTest)

waLBerla_add_test_executable( P23DSmoothConvergenceTest P23DSmoothConvergenceTest.cpp )
target_link_libraries       ( P23DSmoothConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P23DSmoothConvergenceTest)

waLBerla_add_test_executable( P2JacobiConvergenceTest P2JacobiConvergenceTest.cpp )
target_link_libraries       ( P2JacobiConvergenceTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2JacobiConvergenceTest)

waLBerla_add_test_executable( P2MinResConvergenceTest P2MinResConvergenceTest.cpp )
target_link_libraries       ( P2MinResConvergenceTest hyteg walberla::core constant_stencil_operator )
walberla_execute_test(NAME P2MinResConvergenceTest)
waLBerla_execute_test(NAME P2MinResConvergenceTestMPI COMMAND $<TARGET_FILE:P2MinResConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P2GMGHRefinementConvergenceTest P2GMGHRefinementConvergenceTest.cpp )
target_link_libraries       ( P2GMGHRefinementConvergenceTest hyteg walberla::core constant_stencil_operator )
walberla_execute_test(NAME P2GMGHRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGHRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGHRefinementConvergenceTest> PROCESSES 2)

waLBerla_add_test_executable( P2GMGPRefinementConvergenceTest P2GMGPRefinementConvergenceTest.cpp )
target_link_libraries       ( P2GMGPRefinementConvergenceTest hyteg walberla::core constant_stencil_operator )
walberla_execute_test(NAME P2GMGPRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGPRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGPRefinementConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_add_test_executable( P2GSConvergenceTest P2GSConvergenceTest.cpp )
    target_link_libraries       ( P2GSConvergenceTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P2GSConvergenceTest)
endif()

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_add_test_executable( P1ChebyshevSmootherConvergenceTest P1ChebyshevSmootherConvergenceTest.cpp )
    target_link_libraries       ( P1ChebyshevSmootherConvergenceTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME P1ChebyshevSmootherConvergenceTest)
endif()

waLBerla_add_test_executable( P2UnsteadyDiffusion2DTest P2UnsteadyDiffusion2DTest.cpp )
target_link_libraries       ( P2UnsteadyDiffusion2DTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME P2UnsteadyDiffusion2DTest)

waLBerla_add_test_executable( P2UnsteadyDiffusionBlending2DTest P2UnsteadyDiffusionBlending2DTest.cpp )
target_link_libraries       ( P2UnsteadyDiffusionBlending2DTest hyteg walberla::core constant_stencil_operator )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME P2UnsteadyDiffusionBlending2DTest)
    endif ()
endif ()

