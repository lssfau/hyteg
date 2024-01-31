if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/2019_constanta")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME MultigridStudies
        FILES MultigridStudies.cpp
        DEPENDS hyteg core sqlite mixedOperator)

add_subdirectory(2019_constanta)
add_subdirectory(2019_tme)
add_subdirectory(2019_supermuc)
add_subdirectory(2020_agglomeration)
add_subdirectory(2020_supermuc)