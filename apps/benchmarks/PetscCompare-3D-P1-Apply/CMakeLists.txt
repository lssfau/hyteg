if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )
waLBerla_link_files_to_builddir( *.py )

waLBerla_add_executable( NAME PetscCompare-3D-P1-Apply
        FILES PetscCompare-3D-P1-Apply.cpp
        DEPENDS hyteg core constantStencilOperator)