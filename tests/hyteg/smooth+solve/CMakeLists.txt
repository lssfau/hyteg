## Basic Smoothing and Solving
waLBerla_add_test_executable( SmoothableTest SmoothableTest.cpp )
target_link_libraries       ( SmoothableTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME SmoothableTest)

waLBerla_add_test_executable( VectorToVectorOperatorCGTest VectorToVectorOperatorCGTest.cpp )
target_link_libraries       ( VectorToVectorOperatorCGTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME VectorToVectorOperatorCGTest)

waLBerla_add_test_executable( VectorToVectorOperatorChebyshevTest VectorToVectorOperatorChebyshevTest.cpp )
target_link_libraries       ( VectorToVectorOperatorChebyshevTest hyteg walberla::core mixed_operator )
waLBerla_execute_test(NAME VectorToVectorOperatorChebyshevTest)

