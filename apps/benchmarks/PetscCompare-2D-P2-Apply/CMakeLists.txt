if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME PetscCompare-2D-P2-Apply
        FILES PetscCompare-2D-P2-Apply.cpp
        DEPENDS hyteg core constant_stencil_operator)