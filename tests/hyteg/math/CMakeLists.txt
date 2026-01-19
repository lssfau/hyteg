waLBerla_add_test_executable( IntersectionTest IntersectionTest.cpp )
target_link_libraries       ( IntersectionTest hyteg walberla::core )
waLBerla_execute_test(NAME IntersectionTest)

