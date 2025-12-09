if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

add_executable       ( PetscCompare-3D-P1-Apply PetscCompare-3D-P1-Apply.cpp )
target_link_libraries( PetscCompare-3D-P1-Apply hyteg walberla::core constant_stencil_operator )