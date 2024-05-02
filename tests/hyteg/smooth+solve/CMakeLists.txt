## Basic Smoothing and Solving
waLBerla_compile_test(FILES SmoothableTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME SmoothableTest)

waLBerla_compile_test(FILES VectorToVectorOperatorCGTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorToVectorOperatorCGTest)

waLBerla_compile_test(FILES VectorToVectorOperatorChebyshevTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VectorToVectorOperatorChebyshevTest)

