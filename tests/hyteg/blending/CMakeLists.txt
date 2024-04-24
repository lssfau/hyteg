waLBerla_compile_test(FILES AffineMap3DTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME AffineMap3DTest1 COMMAND $<TARGET_FILE:AffineMap3DTest>)
endif ()

waLBerla_compile_test(FILES AnnulusMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME AnnulusMapTest1 COMMAND $<TARGET_FILE:AnnulusMapTest>)

waLBerla_compile_test(FILES BasicBlendingMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BasicBlendingMapTest)

waLBerla_compile_test(FILES IcosahedralShellMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME IcosahedralShellMapTest COMMAND $<TARGET_FILE:IcosahedralShellMapTest>)

waLBerla_compile_test(FILES InverseBlending+PrimitiveIDFindingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME InverseBlending+PrimitiveIDFindingTest)

waLBerla_compile_test(FILES PolygonDistanceTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME PolygonDistanceTest)

waLBerla_compile_test(FILES TokamakMapTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME TokamakMapTest)

if (NOT WALBERLA_CXX_COMPILER_IS_INTEL)
    waLBerla_compile_test(FILES P2P1StokesElementwiseAnnulusBlendingTest.cpp DEPENDS hyteg core mixed_operator)
    waLBerla_execute_test(NAME P2P1StokesElementwiseAnnulusBlendingTest COMMAND $<TARGET_FILE:P2P1StokesElementwiseAnnulusBlendingTest>)
endif ()

waLBerla_compile_test(FILES EvalDFinvDFTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EvalDFinvDFTest)

