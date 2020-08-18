waLBerla_compile_test(FILES MMOC2DCircularConvectionTest.cpp DEPENDS hyteg core convection_particles )
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DCircularConvectionTest)
endif()

waLBerla_compile_test(FILES MMOC2DCircularConvectionBlendingTest.cpp DEPENDS hyteg core convection_particles )
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DCircularConvectionBlendingTest)
endif()

waLBerla_compile_test(FILES MMOC3DCircularConvectionTest.cpp DEPENDS hyteg core convection_particles )
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC3DCircularConvectionTest LABELS longrun)
endif()

waLBerla_compile_test(FILES MMOC2DSwirlingFlowConvectionTest.cpp DEPENDS hyteg core convection_particles )
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DSwirlingFlowConvectionTest)
endif()

waLBerla_compile_test(FILES ConvectionParticlesCouplingTest.cpp DEPENDS core convection_particles )
waLBerla_execute_test(NAME ConvectionParticlesCouplingTest)

waLBerla_compile_test(FILES P2UnsteadyConvectionDiffusion2DTest.cpp DEPENDS hyteg core convection_particles )
waLBerla_execute_test(NAME P2UnsteadyConvectionDiffusion2DTest)