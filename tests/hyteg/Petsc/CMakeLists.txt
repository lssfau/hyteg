if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_add_test_executable( PetscCreateFromLinCombTest PetscCreateFromLinCombTest.cpp )
    target_link_libraries       ( PetscCreateFromLinCombTest hyteg walberla::core )
    waLBerla_execute_test(NAME PetscCreateFromLinCombTest COMMAND $<TARGET_FILE:PetscCreateFromLinCombTest>)

    waLBerla_add_test_executable( PetscMatrixAssemblyTest PetscMatrixAssemblyTest.cpp )
    target_link_libraries       ( PetscMatrixAssemblyTest hyteg walberla::core mixed_operator opgen-diffusion )
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest>)
    waLBerla_execute_test(NAME PetscMatrixAssemblyTest2 COMMAND $<TARGET_FILE:PetscMatrixAssemblyTest> PROCESSES 2)

    waLBerla_add_test_executable( PetscTest PetscTest.cpp )
    target_link_libraries       ( PetscTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME PetscTest COMMAND $<TARGET_FILE:PetscTest>)
    waLBerla_execute_test(NAME PetscTest2 COMMAND $<TARGET_FILE:PetscTest> PROCESSES 2)

    waLBerla_add_test_executable( PetscApplyDirichletBoundaryUIDTest PetscApplyDirichletBoundaryUIDTest.cpp )
    target_link_libraries       ( PetscApplyDirichletBoundaryUIDTest hyteg walberla::core constant_stencil_operator )
    waLBerla_execute_test(NAME PetscApplyDirichletBoundaryUIDTest)
endif ()

if (HYTEG_BUILD_WITH_PETSC AND HYTEG_PETSC_BUILT_WITH_HDF5)
    waLBerla_add_test_executable( PetscHDF5SaveLoadTest PetscHDF5SaveLoadTest.cpp )
    target_link_libraries       ( PetscHDF5SaveLoadTest hyteg walberla::core )
    waLBerla_execute_test(NAME PetscHDF5SaveLoadTest COMMAND $<TARGET_FILE:PetscHDF5SaveLoadTest>)
endif ()

