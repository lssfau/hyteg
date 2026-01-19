waLBerla_add_test_executable( EdgeDoFAdditiveCommunicationTest EdgeDoFAdditiveCommunicationTest.cpp )
target_link_libraries       ( EdgeDoFAdditiveCommunicationTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFAdditiveCommunicationTest)

waLBerla_add_test_executable( EdgeDoFMemoryTest EdgeDoFMemoryTest.cpp )
target_link_libraries       ( EdgeDoFMemoryTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFMemoryTest)

waLBerla_add_test_executable( EdgeDoFFunctionTest EdgeDoFFunctionTest.cpp )
target_link_libraries       ( EdgeDoFFunctionTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFFunctionTest)

waLBerla_add_test_executable( EdgeDoFEnumerateTest EdgeDoFEnumerateTest.cpp )
target_link_libraries       ( EdgeDoFEnumerateTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFEnumerateTest)

waLBerla_add_test_executable( EdgeDoFCommTest EdgeDoFCommTest.cpp )
target_link_libraries       ( EdgeDoFCommTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFCommTest)

waLBerla_add_test_executable( EdgeDoFInterpolationTest EdgeDoFInterpolationTest.cpp )
target_link_libraries       ( EdgeDoFInterpolationTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFInterpolationTest)

waLBerla_add_test_executable( EdgeDoFFunction3DTest EdgeDoFFunction3DTest.cpp )
target_link_libraries       ( EdgeDoFFunction3DTest hyteg walberla::core )
waLBerla_execute_test(NAME EdgeDoFFunction3DTest)

