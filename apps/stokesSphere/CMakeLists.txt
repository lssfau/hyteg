if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_add_executable( NAME StokesSphere
        FILES  StokesSphere.cpp
        DEPENDS tinyhhg_core)