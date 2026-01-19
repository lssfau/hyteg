waLBerla_add_test_executable( SphereToolsTest SphereToolsTest.cpp )
target_link_libraries       ( SphereToolsTest hyteg walberla::core )
waLBerla_execute_test(NAME SphereToolsTest1 COMMAND $<TARGET_FILE:SphereToolsTest>)
waLBerla_execute_test(NAME SphereToolsTest2 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 2)
waLBerla_execute_test(NAME SphereToolsTest4 COMMAND $<TARGET_FILE:SphereToolsTest> PROCESSES 4)

waLBerla_add_test_executable( closestPointTriangleTetrahedron3DTest closestPointTriangleTetrahedron3DTest.cpp )
target_link_libraries       ( closestPointTriangleTetrahedron3DTest hyteg walberla::core )
waLBerla_execute_test(NAME closestPointTriangleTetrahedron3DTest)