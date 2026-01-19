## adaptive refinement w/ hanging nodes
waLBerla_add_test_executable( ARHangingNodesStorageTest ARHangingNodesStorageTest.cpp )
target_link_libraries       ( ARHangingNodesStorageTest hyteg walberla::core )
waLBerla_execute_test(NAME ARHangingNodesStorageTest)

waLBerla_add_test_executable( MacroVolumeRefinementMicroIdxConversionTest MacroVolumeRefinementMicroIdxConversionTest.cpp )
target_link_libraries       ( MacroVolumeRefinementMicroIdxConversionTest hyteg walberla::core )
waLBerla_execute_test(NAME MacroVolumeRefinementMicroIdxConversionTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( P1AdaptiveRefinementEqualityTest P1AdaptiveRefinementEqualityTest.cpp )
    target_link_libraries       ( P1AdaptiveRefinementEqualityTest hyteg constant_stencil_operator )
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME P1AdaptiveRefinementEqualityTest)
    endif()

    waLBerla_add_test_executable( DGAdaptiveRefinementEqualityTest DGAdaptiveRefinementEqualityTest.cpp )
    target_link_libraries       ( DGAdaptiveRefinementEqualityTest hyteg constant_stencil_operator )
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME DGAdaptiveRefinementEqualityTest)
    endif()
endif ()
