waLBerla_link_files_to_builddir( *.prm )

waLBerla_compile_test(FILES ElementwiseDivKGradCGConvergenceTest.cpp DEPENDS hyteg core opgen-div_k_grad)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradCGConvergenceTest)
endif ()

waLBerla_compile_test(FILES ElementwiseDivKGradNeumannCGConvergenceTest.cpp DEPENDS hyteg core opgen-div_k_grad)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME ElementwiseDivKGradNeumannCGConvergenceTest)
endif ()

waLBerla_compile_test(FILES P2P1StokesPreconditionersConvergenceTest.cpp DEPENDS hyteg core opgen-k_mass opgen-composites-stokes constant_stencil_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME P2P1StokesPreconditionersConvergenceTest)
    waLBerla_execute_test(NAME P2P1StokesPreconditionersConvergenceTestLongrun COMMAND $<TARGET_FILE:P2P1StokesPreconditionersConvergenceTest> --longrun LABELS longrun)
endif ()

waLBerla_compile_test(FILES P1VectorLaplaceCGTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P1VectorLaplaceCGTest)

waLBerla_compile_test(FILES P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P1StokesMinResConvergenceTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1P1StokesUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1P1StokesUzawaConvergenceTest)
endif ()

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_compile_test(FILES P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P1Stokes3DMinResConvergenceTest)

    waLBerla_compile_test(FILES P2P1Stokes3DMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes3DMinResConvergenceTest)

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES P1P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
        waLBerla_execute_test(NAME P1P1Stokes3DUzawaConvergenceTest)
    endif ()

    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_compile_test(FILES P2P1Stokes3DUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
        waLBerla_execute_test(NAME P2P1Stokes3DUzawaConvergenceTest)
    endif ()
endif ()

waLBerla_compile_test(FILES P2ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P2ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P2ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES ElementwiseEpsilonMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator opgen-composites-stokes)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTest)
        waLBerla_execute_test(NAME ElementwiseEpsilonMinResConvergenceTestLongrun COMMAND $<TARGET_FILE:ElementwiseEpsilonMinResConvergenceTest> --longrun LABELS longrun)
    endif ()
endif ()

waLBerla_compile_test(FILES P2P1StokesMinResConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1StokesMinResConvergenceTest)

waLBerla_compile_test(FILES P2P1StokesSchurCGConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1StokesSchurCGConvergenceTest)

waLBerla_compile_test(FILES P2P1ElementwiseUzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1ElementwiseUzawaConvergenceTest)

waLBerla_compile_test(FILES P2P1UzawaConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME P2P1UzawaConvergenceTest)

waLBerla_compile_test(FILES StokesVarViscPopovTest.cpp DEPENDS hyteg core mixed_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME StokesVarViscPopovTest)
    endif ()
endif ()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES P2P1Stokes2DGKBConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1Stokes2DGKBConvergenceTest)
endif()

if( HYTEG_BUILD_WITH_PETSC )
    waLBerla_compile_test(FILES P2P1ChannelGKBConvergenceTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1ChannelGKBConvergenceTest)
endif()

waLBerla_compile_test(FILES P1CGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1CGConvergenceTest)
waLBerla_execute_test(NAME P1CGConvergenceTestMPI COMMAND $<TARGET_FILE:P1CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1CG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1CG3DConvergenceTest)
waLBerla_execute_test(NAME P1CG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1CG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1GMGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMGConvergenceTest)
waLBerla_execute_test(NAME P1GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1FASConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1FASConvergenceTest)
waLBerla_execute_test(NAME P1FASConvergenceTestMPI COMMAND $<TARGET_FILE:P1FASConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1TauExtrapolationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTest)
waLBerla_execute_test(NAME P1TauExtrapolationConvergenceTestMPI COMMAND $<TARGET_FILE:P1TauExtrapolationConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1GMG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMG3DConvergenceTest)
waLBerla_execute_test(NAME P1GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P1GMG3DConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P2GMGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMGConvergenceTest)
waLBerla_execute_test(NAME P2GMGConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P2GMG3DConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMG3DConvergenceTest)
waLBerla_execute_test(NAME P2GMG3DConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMG3DConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES P1GSConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1GSConvergenceTest)
endif()

waLBerla_compile_test(FILES P1JacobiConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1JacobiConvergenceTest)

waLBerla_compile_test(FILES P1MinResConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1MinResConvergenceTest)

waLBerla_compile_test(FILES P2CGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2CGConvergenceTest)
waLBerla_execute_test(NAME P2CGConvergenceTestMPI COMMAND $<TARGET_FILE:P2CGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P1ElementwiseCGConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTest)
waLBerla_execute_test(NAME P1ElementwiseCGConvergenceTestMPI COMMAND $<TARGET_FILE:P1ElementwiseCGConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P2P1FGMRESStokesTest.cpp DEPENDS hyteg core constant_stencil_operator mixed_operator)
waLBerla_execute_test(NAME P2P1FGMRESStokesTest)

waLBerla_compile_test(FILES P2P1GMRESStokesTest.cpp DEPENDS hyteg core constant_stencil_operator mixed_operator)
waLBerla_execute_test(NAME P2P1GMRESStokesTest)

waLBerla_compile_test(FILES P2GMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2GMRESConvergenceTest)

waLBerla_compile_test(FILES P2FGMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2FGMRESConvergenceTest)

waLBerla_compile_test(FILES P1GMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1GMRESConvergenceTest)

waLBerla_compile_test(FILES P1FGMRESConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P1FGMRESConvergenceTest)

waLBerla_compile_test(FILES P2SmoothConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2SmoothConvergenceTest)

waLBerla_compile_test(FILES P23DSmoothConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P23DSmoothConvergenceTest)

waLBerla_compile_test(FILES P2JacobiConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2JacobiConvergenceTest)

waLBerla_compile_test(FILES P2MinResConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2MinResConvergenceTest)
waLBerla_execute_test(NAME P2MinResConvergenceTestMPI COMMAND $<TARGET_FILE:P2MinResConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P2GMGHRefinementConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2GMGHRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGHRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGHRefinementConvergenceTest> PROCESSES 2)

waLBerla_compile_test(FILES P2GMGPRefinementConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
walberla_execute_test(NAME P2GMGPRefinementConvergenceTest)
waLBerla_execute_test(NAME P2GMGPRefinementConvergenceTestMPI COMMAND $<TARGET_FILE:P2GMGPRefinementConvergenceTest> PROCESSES 2)

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES P2GSConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P2GSConvergenceTest)
endif()

if(WALBERLA_DOUBLE_ACCURACY)
    waLBerla_compile_test(FILES P1ChebyshevSmootherConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME P1ChebyshevSmootherConvergenceTest)
endif()

waLBerla_compile_test(FILES P2UnsteadyDiffusion2DTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME P2UnsteadyDiffusion2DTest)

waLBerla_compile_test(FILES P2UnsteadyDiffusionBlending2DTest.cpp DEPENDS hyteg core constant_stencil_operator)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (HYTEG_BUILD_WITH_PETSC)
        waLBerla_execute_test(NAME P2UnsteadyDiffusionBlending2DTest)
    endif ()
endif ()

