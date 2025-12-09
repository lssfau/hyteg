if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( VectorFunctionPetscTest VectorFunctionPetscTest.cpp )
    target_link_libraries       ( VectorFunctionPetscTest hyteg walberla::core )
    waLBerla_execute_test(NAME VectorFunctionPetscTest)
endif ()

waLBerla_add_test_executable( VectorFunctionBasicTest VectorFunctionBasicTest.cpp )
target_link_libraries       ( VectorFunctionBasicTest hyteg walberla::core constant_stencil_operator )
waLBerla_execute_test(NAME VectorFunctionBasicTest)

