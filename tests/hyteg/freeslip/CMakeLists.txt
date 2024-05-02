waLBerla_compile_test(FILES ProjectNormalTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ProjectNormalTest)
waLBerla_execute_test(NAME ProjectNormalTestMPI COMMAND $<TARGET_FILE:ProjectNormalTest> PROCESSES 2)

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES FreeslipPetscApplyTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME FreeslipPetscApplyTest)
endif ()

waLBerla_compile_test(FILES FreeslipRectangularChannelTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME FreeslipRectangularChannelTest)
waLBerla_execute_test(NAME FreeslipRectangularChannelTestMPI COMMAND $<TARGET_FILE:FreeslipRectangularChannelTest> PROCESSES 2)

