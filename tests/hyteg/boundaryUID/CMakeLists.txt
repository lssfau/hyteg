waLBerla_add_test_executable( BoundaryUIDTest BoundaryUIDTest.cpp )
target_link_libraries       ( BoundaryUIDTest hyteg walberla::core )
waLBerla_execute_test(NAME BoundaryUIDTest)
waLBerla_add_test_executable( HollowBodyUIDTest HollowBodyUIDTest.cpp )
target_link_libraries       ( HollowBodyUIDTest hyteg walberla::core )
waLBerla_execute_test(NAME HollowBodyUIDTest)

