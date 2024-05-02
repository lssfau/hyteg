waLBerla_compile_test(FILES EdgeDoFAdditiveCommunicationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFAdditiveCommunicationTest)

waLBerla_compile_test(FILES EdgeDoFMemoryTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFMemoryTest)

waLBerla_compile_test(FILES EdgeDoFFunctionTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFFunctionTest)

waLBerla_compile_test(FILES EdgeDoFEnumerateTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFEnumerateTest)

waLBerla_compile_test(FILES EdgeDoFCommTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFCommTest)

waLBerla_compile_test(FILES EdgeDoFInterpolationTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFInterpolationTest)

waLBerla_compile_test(FILES EdgeDoFFunction3DTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME EdgeDoFFunction3DTest)

