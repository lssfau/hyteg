## adaptive refinement w/ hanging nodes
waLBerla_compile_test(FILES ARHangingNodesStorageTest.cpp)
waLBerla_execute_test(NAME ARHangingNodesStorageTest)

waLBerla_compile_test(FILES MacroVolumeRefinementMicroIdxConversionTest.cpp)
waLBerla_execute_test(NAME MacroVolumeRefinementMicroIdxConversionTest)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES P1AdaptiveRefinementEqualityTest.cpp DEPENDS hyteg constant_stencil_operator)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME P1AdaptiveRefinementEqualityTest)
    endif()

    waLBerla_compile_test(FILES DGAdaptiveRefinementEqualityTest.cpp DEPENDS hyteg constant_stencil_operator)
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        waLBerla_execute_test(NAME DGAdaptiveRefinementEqualityTest)
    endif()
endif ()

