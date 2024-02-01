if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()



waLBerla_link_files_to_builddir( *.prm )

waLBerla_link_files_to_builddir( *.sh )

waLBerla_add_executable( NAME 11_CouetteFlow2D
        FILES  11_CouetteFlow2D.cpp
        DEPENDS hyteg mixed_operator)
