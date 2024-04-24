waLBerla_compile_test(FILES BoundaryUIDTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME BoundaryUIDTest)
waLBerla_compile_test(FILES HollowBodyUIDTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME HollowBodyUIDTest)

