waLBerla_compile_test(FILES json/BasicJSONReadTest.cpp DEPENDS hyteg core )
waLBerla_execute_test(NAME BasicJSONReadTest)

waLBerla_compile_test(FILES plates/PlateVelocityComputationTest.cpp DEPENDS terraneo hyteg core )
waLBerla_execute_test(NAME PlateVelocityComputationTest)

waLBerla_compile_test(FILES initialisation/InitialisationTest.cpp DEPENDS terraneo hyteg core )
waLBerla_execute_test(NAME InitialisationTest)

waLBerla_compile_test(FILES DataIO/TerraNeoDataIOTest.cpp DEPENDS terraneo hyteg core )
waLBerla_execute_test(NAME TerraNeoDataIOTest)
waLBerla_link_files_to_builddir(*.prm)