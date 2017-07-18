if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/meshes)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/meshes")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/param)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/param")
endif()

if(NOT EXISTS ${CMAKE_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir(../data/meshes/*.msh)
waLBerla_link_files_to_builddir(../data/param/*.prm)

waLBerla_add_executable( NAME tinytest
        FILES test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_cg
        FILES cg_test.cpp
        DEPENDS tinyhhg_core fmt core)

waLBerla_add_executable( NAME tinytest_cg_bubble
        FILES cg_bubble_test.cpp
        DEPENDS tinyhhg_core fmt core)

waLBerla_add_executable( NAME tinytest_chorin
        FILES chorin_test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_minres
        FILES minres_test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_block_cg
        FILES block_cg_test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_stokes_stab
        FILES stokes_stab.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_fullmg
        FILES fullmg_test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME cg_P1_new
        FILES cg_P1_new.cpp
        DEPENDS tinyhhg_core fmt)

#add_executable (tinytest_fmg fmg_test.cpp)
#target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})
#
#add_executable (tinytest_chorin chorin_test.cpp)
#target_link_libraries (tinytest_chorin tinyhhg fmt ${MPI_CXX_LIBRARIES})
