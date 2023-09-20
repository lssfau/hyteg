waLBerla_compile_test(FILES UnresolvedParticlesCouplingTest.cpp DEPENDS core unresolved_particles coupling_hyteg_unresolved_particles )
waLBerla_execute_test(NAME UnresolvedParticlesCouplingTest)

waLBerla_compile_test(FILES UnresolvedParticlesSimpleTransportTest.cpp DEPENDS core unresolved_particles coupling_hyteg_unresolved_particles )
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME UnresolvedParticlesSimpleTransportTest)
    waLBerla_execute_test(NAME UnresolvedParticlesSimpleTransportTest2 COMMAND $<TARGET_FILE:UnresolvedParticlesSimpleTransportTest> PROCESSES 2 )
    waLBerla_execute_test(NAME UnresolvedParticlesSimpleTransportTest4 COMMAND $<TARGET_FILE:UnresolvedParticlesSimpleTransportTest> PROCESSES 4 )
endif()