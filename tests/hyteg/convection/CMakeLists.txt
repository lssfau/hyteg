waLBerla_add_test_executable( Upwind2DCircularConvectionTest Upwind2DCircularConvectionTest.cpp )
target_link_libraries       ( Upwind2DCircularConvectionTest hyteg walberla::core )
waLBerla_add_test_executable( AlgebraicUpwind3DCircularConvectionTest AlgebraicUpwind3DCircularConvectionTest.cpp )
target_link_libraries       ( AlgebraicUpwind3DCircularConvectionTest hyteg walberla::core constant_stencil_operator )

