if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

add_executable       ( StokesFreeSlip StokesFreeSlip.cpp )
target_link_libraries( StokesFreeSlip hyteg mixed_operator )

add_executable       ( ProjectNormalDemo ProjectNormalDemo.cpp )
target_link_libraries( ProjectNormalDemo hyteg )
