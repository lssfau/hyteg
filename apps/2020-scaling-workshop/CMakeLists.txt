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
add_executable       ( Scaling_Workshop_01_Cube Scaling_Workshop_01_Cube.cpp Helpers.cpp )
target_link_libraries( Scaling_Workshop_01_Cube hyteg walberla::sqlite mixed_operator )

add_executable       ( Scaling_Workshop_02_Spherical_Shell Scaling_Workshop_02_Spherical_Shell.cpp Helpers.cpp )
target_link_libraries( Scaling_Workshop_02_Spherical_Shell hyteg walberla::sqlite mixed_operator )
endif()

add_subdirectory(scriptGeneration)