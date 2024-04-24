if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES VectorFunctionPetscTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME VectorFunctionPetscTest)
endif ()

waLBerla_compile_test(FILES VectorFunctionBasicTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME VectorFunctionBasicTest)

