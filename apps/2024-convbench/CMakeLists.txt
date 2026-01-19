if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

# include(CMakePrintHelpers)
# cmake_print_variables(Python3_VERSION)
# cmake_print_variables(Python3_INCLUDE_DIRS)

# include_directories(${Python3_INCLUDE_DIRS})

waLBerla_link_files_to_builddir( *.py )
waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *sh )

add_executable       ( 2DStokesFSConvTest 2DStokesFSConvTest.cpp )
target_link_libraries( 2DStokesFSConvTest hyteg mixed_operator )

add_executable       ( 3DStokesFSConvTest 3DStokesFSConvTest.cpp )
target_link_libraries( 3DStokesFSConvTest hyteg mixed_operator )
