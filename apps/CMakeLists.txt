if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/meshes)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/meshes")
endif()

waLBerla_link_files_to_builddir(../data/meshes/*.msh)

waLBerla_add_executable( NAME tinytest
        FILES test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_cg
        FILES cg_test.cpp
        DEPENDS tinyhhg_core fmt core)

waLBerla_add_executable( NAME tinytest_fmg
        FILES fmg_test.cpp
        DEPENDS tinyhhg_core fmt)

waLBerla_add_executable( NAME tinytest_chorin
        FILES chorin_test.cpp
        DEPENDS tinyhhg_core fmt)

#add_executable (tinytest_fmg fmg_test.cpp)
#target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})
#
#add_executable (tinytest_chorin chorin_test.cpp)
#target_link_libraries (tinytest_chorin tinyhhg fmt ${MPI_CXX_LIBRARIES})
