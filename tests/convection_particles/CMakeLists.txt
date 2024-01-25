waLBerla_compile_test(FILES MMOC2DCircularConvectionTest.cpp DEPENDS hyteg core convection_particles constant_stencil_operator)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DCircularConvectionTest)
    waLBerla_execute_test(NAME MMOC2DCircularConvectionTestMPI2 COMMAND $<TARGET_FILE:MMOC2DCircularConvectionTest> PROCESSES 2 )
endif()

waLBerla_compile_test(FILES MMOC2DCircularConvectionBlendingTest.cpp DEPENDS hyteg core convection_particles constant_stencil_operator)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DCircularConvectionBlendingTest)
    waLBerla_execute_test(NAME MMOC2DCircularConvectionBlendingTestMPI2 COMMAND $<TARGET_FILE:MMOC2DCircularConvectionBlendingTest> PROCESSES 2 )
    waLBerla_execute_test(NAME MMOC2DCircularConvectionBlendingTestMPI4 COMMAND $<TARGET_FILE:MMOC2DCircularConvectionBlendingTest> PROCESSES 4 )
endif()

waLBerla_compile_test(FILES MMOC3DCircularConvectionTest.cpp DEPENDS hyteg core convection_particles constant_stencil_operator)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC3DCircularConvectionTest LABELS longrun)
    waLBerla_execute_test(NAME MMOC3DCircularConvectionTestMPI2 COMMAND $<TARGET_FILE:MMOC3DCircularConvectionTest> PROCESSES 2 )
    waLBerla_execute_test(NAME MMOC3DCircularConvectionTestMPI4 COMMAND $<TARGET_FILE:MMOC3DCircularConvectionTest> PROCESSES 4 )
endif()

waLBerla_compile_test(FILES MMOC2DSwirlingFlowConvectionTest.cpp DEPENDS hyteg core convection_particles constant_stencil_operator)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test(NAME MMOC2DSwirlingFlowConvectionTest)
    waLBerla_execute_test(NAME MMOC2DSwirlingFlowConvectionTestMPI2 COMMAND $<TARGET_FILE:MMOC2DSwirlingFlowConvectionTest> PROCESSES 2 )

endif()

waLBerla_compile_test(FILES ConvectionParticlesCouplingTest.cpp DEPENDS core convection_particles )
waLBerla_execute_test(NAME ConvectionParticlesCouplingTest)

waLBerla_compile_test(FILES P2UnsteadyConvectionDiffusion2DTest.cpp DEPENDS hyteg core convection_particles constant_stencil_operator)
waLBerla_execute_test(NAME P2UnsteadyConvectionDiffusion2DTest)
waLBerla_execute_test(NAME P2UnsteadyConvectionDiffusion2DTestMPI2 COMMAND $<TARGET_FILE:P2UnsteadyConvectionDiffusion2DTest> PROCESSES 2 )
waLBerla_execute_test(NAME P2UnsteadyConvectionDiffusion2DTestMPI4 COMMAND $<TARGET_FILE:P2UnsteadyConvectionDiffusion2DTest> PROCESSES 4 )