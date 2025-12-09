if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/2019_constanta")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( MultigridStudies MultigridStudies.cpp )
target_link_libraries( MultigridStudies hyteg walberla::core walberla::sqlite mixed_operator )

add_subdirectory(2019_constanta)
add_subdirectory(2019_tme)
add_subdirectory(2019_supermuc)
add_subdirectory(2020_agglomeration)
add_subdirectory(2020_supermuc)