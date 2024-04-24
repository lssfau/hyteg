if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES PetscCreateFromLinCombTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PetscCreateFromLinCombTest COMMAND $<TARGET_FILE:PetscCreateFromLinCombTest>)

    waLBerla_compile_test(FILES PetscMatrixAssemblyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest>)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest2 COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest> PROCESSES 2)

    waLBerla_compile_test(FILES PetscTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest>)
    waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2)

    waLBerla_compile_test(FILES PetscApplyDirichletBoundaryUIDTest.cpp DEPENDS hyteg core constant_stencil_operator)
    waLBerla_execute_test(NAME PetscApplyDirichletBoundaryUIDTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC AND HYTEG_PETSC_BUILT_WITH_HDF5)
    waLBerla_compile_test(FILES PetscHDF5SaveLoadTest.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME PetscHDF5SaveLoadTest COMMAND $<TARGET_FILE:PetscHDF5SaveLoadTest>)
endif ()

