waLBerla_add_test_executable( BasicJSONReadTest json/BasicJSONReadTest.cpp )
target_link_libraries       ( BasicJSONReadTest hyteg walberla::core )
waLBerla_execute_test(NAME BasicJSONReadTest)

waLBerla_add_test_executable( PlateVelocityComputationTest plates/PlateVelocityComputationTest.cpp )
target_link_libraries       ( PlateVelocityComputationTest terraneo hyteg walberla::core )
waLBerla_execute_test(NAME PlateVelocityComputationTest)

waLBerla_add_test_executable( InitialisationTest initialisation/InitialisationTest.cpp )
target_link_libraries       ( InitialisationTest terraneo hyteg walberla::core )
waLBerla_execute_test(NAME InitialisationTest)

waLBerla_add_test_executable( TerraNeoParameterTest_v0_1 parameter/TerraNeoParameterTest_v0_1.cpp )
target_link_libraries       ( TerraNeoParameterTest_v0_1 terraneo hyteg walberla::core )
waLBerla_execute_test(NAME TerraNeoParameterTest_v0_1)

waLBerla_add_test_executable( TerraNeoParameterTest_v0_2 parameter/TerraNeoParameterTest_v0_2.cpp )
target_link_libraries       ( TerraNeoParameterTest_v0_2 terraneo hyteg walberla::core )
waLBerla_execute_test(NAME TerraNeoParameterTest_v0_2)

waLBerla_add_test_executable( RadialOutputTest radialoutput/RadialOutputTest.cpp )
target_link_libraries       ( RadialOutputTest terraneo hyteg walberla::core )
waLBerla_execute_test(NAME RadialOutputTest)

waLBerla_add_test_executable( MeshTransferWithParticlesTest utils/MeshTransferWithParticlesTest.cpp )
target_link_libraries       ( MeshTransferWithParticlesTest terraneo hyteg walberla::core convection_particles opgen-mass )
waLBerla_execute_test(NAME MeshTransferWithParticlesTest1 COMMAND $<TARGET_FILE:MeshTransferWithParticlesTest>)
waLBerla_execute_test(NAME MeshTransferWithParticlesTest4 COMMAND $<TARGET_FILE:MeshTransferWithParticlesTest> PROCESSES 4)
waLBerla_execute_test(NAME MeshTransferWithParticlesTest12 COMMAND $<TARGET_FILE:MeshTransferWithParticlesTest> PROCESSES 12)

add_subdirectory( solvers )
