add_executable (tinytest test.cpp)
target_link_libraries (tinytest tinyhhg fmt ${MPI_CXX_LIBRARIES})

add_executable (tinytest_cg cg_test.cpp)
target_link_libraries (tinytest_cg tinyhhg fmt ${MPI_CXX_LIBRARIES})

add_executable (tinytest_fmg fmg_test.cpp)
target_link_libraries (tinytest_fmg tinyhhg fmt ${MPI_CXX_LIBRARIES} ${LIKWID_LIB})

add_executable (tinytest_chorin chorin_test.cpp)
target_link_libraries (tinytest_chorin tinyhhg fmt ${MPI_CXX_LIBRARIES})
