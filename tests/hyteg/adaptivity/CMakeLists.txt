waLBerla_compile_test(FILES PrimitiveMigrationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveMigrationTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationTest>)
waLBerla_execute_test(NAME PrimitiveMigrationTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationTest8 COMMAND $<TARGET_FILE:PrimitiveMigrationTest> PROCESSES 8)

waLBerla_compile_test(FILES ParallelRoundRobinTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME ParallelRoundRobinTest1 COMMAND $<TARGET_FILE:ParallelRoundRobinTest>)
waLBerla_execute_test(NAME ParallelRoundRobinTest3 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 3)
waLBerla_execute_test(NAME ParallelRoundRobinTest8 COMMAND $<TARGET_FILE:ParallelRoundRobinTest> PROCESSES 8)

waLBerla_compile_test(FILES PrimitiveStorageCopyTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest1 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest>)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest3 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveStorageCopyTest8 COMMAND $<TARGET_FILE:PrimitiveStorageCopyTest> PROCESSES 8)

if (WALBERLA_BUILD_WITH_PARMETIS)
    waLBerla_compile_test(FILES ParMetisTestHyteg.cpp DEPENDS hyteg core)
    waLBerla_execute_test(NAME ParMetisTestHyteg1 COMMAND $<TARGET_FILE:ParMetisTestHyteg>)
    waLBerla_execute_test(NAME ParMetisTestHyteg3 COMMAND $<TARGET_FILE:ParMetisTestHyteg> PROCESSES 3)
    waLBerla_execute_test(NAME ParMetisTestHyteg8 COMMAND $<TARGET_FILE:ParMetisTestHyteg> PROCESSES 8)
endif ()

waLBerla_compile_test(FILES PrimitiveMigrationMatMulTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest>)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationMatMulTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationMatMulTest> PROCESSES 4)

waLBerla_compile_test(FILES PrimitiveMigrationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest1 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest>)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest3 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 3)
waLBerla_execute_test(NAME PrimitiveMigrationConvergenceTest4 COMMAND $<TARGET_FILE:PrimitiveMigrationConvergenceTest> PROCESSES 4)

waLBerla_compile_test(FILES AgglomerationMatMulTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME AgglomerationMatMulTest1 COMMAND $<TARGET_FILE:AgglomerationMatMulTest>)
waLBerla_execute_test(NAME AgglomerationMatMulTest2 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationMatMulTest3 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationMatMulTest4 COMMAND $<TARGET_FILE:AgglomerationMatMulTest> PROCESSES 4 )

waLBerla_compile_test(FILES AgglomerationConvergenceTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME AgglomerationConvergenceTest1 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest>)
waLBerla_execute_test(NAME AgglomerationConvergenceTest2 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 2)
waLBerla_execute_test(NAME AgglomerationConvergenceTest3 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 3)
#waLBerla_execute_test(NAME AgglomerationConvergenceTest4 COMMAND $<TARGET_FILE:AgglomerationConvergenceTest> PROCESSES 4 )

if (HYTEG_BUILD_WITH_PETSC)
    waLBerla_compile_test(FILES AgglomerationSupermanFactorizationTest.cpp DEPENDS hyteg core constant_stencil_operator)
    #waLBerla_execute_test(NAME AgglomerationSupermanFactorizationTest4 COMMAND $<TARGET_FILE:AgglomerationSupermanFactorizationTest> PROCESSES 4 )
endif ()

