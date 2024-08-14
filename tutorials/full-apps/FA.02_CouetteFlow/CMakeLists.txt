if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()



waLBerla_link_files_to_builddir( *.prm )

waLBerla_link_files_to_builddir( *.sh )

waLBerla_add_executable( NAME FA.02_CouetteFlow2D
        FILES  FA.02_CouetteFlow2D.cpp
        DEPENDS hyteg mixed_operator)

if(CMAKE_BUILD_TYPE MATCHES "Release")
    waLBerla_execute_test( NO_MODULE_LABEL NAME FA.02_CouetteFlow2D COMMAND $<TARGET_FILE:FA.02_CouetteFlow2D>
            "${CMAKE_BINARY_DIR}//data/param/FA.02_CouetteFlow2D_CI.prm"  )
endif()