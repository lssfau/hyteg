waLBerla_compile_test(FILES RefineCoarseMeshComparisonTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME RefineCoarseMeshComparisonTest)

waLBerla_compile_test(FILES RefineCoarseMeshCGTest.cpp DEPENDS hyteg core constant_stencil_operator)
waLBerla_execute_test(NAME RefineCoarseMeshCGTest)

waLBerla_compile_test(FILES MeshInfoTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MeshInfoTest1 COMMAND $<TARGET_FILE:MeshInfoTest>)
waLBerla_execute_test(NAME MeshInfoTest2 COMMAND $<TARGET_FILE:MeshInfoTest> PROCESSES 2)
waLBerla_compile_test(FILES MeshGenTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME MeshInfoTest3 COMMAND $<TARGET_FILE:MeshGenTest>)
waLBerla_execute_test(NAME MeshInfoTest4 COMMAND $<TARGET_FILE:MeshGenTest> PROCESSES 2)

waLBerla_compile_test(FILES 2DCoarseMeshRefinementTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME 2DCoarseMeshRefinementTest)

waLBerla_compile_test(FILES TokamakLaplaceTest.cpp DEPENDS hyteg core)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    waLBerla_execute_test(NAME TokamakLaplaceTest)
    waLBerla_execute_test(NAME TokamakLaplaceTestLongrun COMMAND $<TARGET_FILE:TokamakLaplaceTest> --longrun LABELS longrun)
endif ()

waLBerla_compile_test(FILES MicroMeshSmokeTest.cpp DEPENDS hyteg core opgen-div_k_grad)
waLBerla_execute_test(NAME MicroMeshSmokeTest)