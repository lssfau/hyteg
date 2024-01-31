if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/db)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/db")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/hawk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hawk")
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/supermuc)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/supermuc")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.sh )

if (HYTEG_BUILD_WITH_PETSC)
waLBerla_add_executable( NAME Scaling_Workshop_01_Cube
        FILES  Scaling_Workshop_01_Cube.cpp Helpers.cpp
        DEPENDS hyteg sqlite mixedOperator)

waLBerla_add_executable( NAME Scaling_Workshop_02_Spherical_Shell
        FILES  Scaling_Workshop_02_Spherical_Shell.cpp Helpers.cpp
        DEPENDS hyteg sqlite mixedOperator)
endif()

add_subdirectory(scriptGeneration)