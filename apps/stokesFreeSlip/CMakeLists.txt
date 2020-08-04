if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME StokesFreeSlip
        FILES StokesFreeSlip.cpp
        DEPENDS hyteg)

waLBerla_add_executable( NAME StokesFreeSlipP2
        FILES StokesFreeSlipP2.cpp
        DEPENDS hyteg)