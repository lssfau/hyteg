waLBerla_add_test_executable( AffineMap3DTest AffineMap3DTest.cpp )
target_link_libraries       ( AffineMap3DTest hyteg walberla::core )
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME AffineMap3DTest1 COMMAND $<TARGET_FILE:AffineMap3DTest>)
endif ()

waLBerla_add_test_executable( AnnulusMapTest AnnulusMapTest.cpp )
target_link_libraries       ( AnnulusMapTest hyteg walberla::core )
waLBerla_execute_test(NAME AnnulusMapTest1 COMMAND $<TARGET_FILE:AnnulusMapTest>)

waLBerla_add_test_executable( BasicBlendingMapTest BasicBlendingMapTest.cpp )
target_link_libraries       ( BasicBlendingMapTest hyteg walberla::core )
waLBerla_execute_test(NAME BasicBlendingMapTest)

waLBerla_add_test_executable( IcosahedralShellMapTest IcosahedralShellMapTest.cpp )
target_link_libraries       ( IcosahedralShellMapTest hyteg walberla::core )
waLBerla_execute_test(NAME IcosahedralShellMapTest COMMAND $<TARGET_FILE:IcosahedralShellMapTest>)

waLBerla_add_test_executable( InverseBlending+PrimitiveIDFindingTest InverseBlending+PrimitiveIDFindingTest.cpp )
target_link_libraries       ( InverseBlending+PrimitiveIDFindingTest hyteg walberla::core )
waLBerla_execute_test(NAME InverseBlending+PrimitiveIDFindingTest)

waLBerla_add_test_executable( PolygonDistanceTest PolygonDistanceTest.cpp )
target_link_libraries       ( PolygonDistanceTest hyteg walberla::core )
waLBerla_execute_test(NAME PolygonDistanceTest)

waLBerla_add_test_executable( TokamakMapTest TokamakMapTest.cpp )
target_link_libraries       ( TokamakMapTest hyteg walberla::core )
waLBerla_execute_test(NAME TokamakMapTest)

if (NOT WALBERLA_CXX_COMPILER_IS_INTEL)
    waLBerla_add_test_executable( P2P1StokesElementwiseAnnulusBlendingTest P2P1StokesElementwiseAnnulusBlendingTest.cpp )
    target_link_libraries       ( P2P1StokesElementwiseAnnulusBlendingTest hyteg walberla::core mixed_operator )
    waLBerla_execute_test(NAME P2P1StokesElementwiseAnnulusBlendingTest COMMAND $<TARGET_FILE:P2P1StokesElementwiseAnnulusBlendingTest>)
endif ()

waLBerla_add_test_executable( EvalDFinvDFTest EvalDFinvDFTest.cpp )
target_link_libraries       ( EvalDFinvDFTest hyteg walberla::core )
waLBerla_execute_test(NAME EvalDFinvDFTest)
