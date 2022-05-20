if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/vtk)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME GMRESApp
        FILES GMRESApp.cpp
        DEPENDS hyteg core )

waLBerla_add_executable( NAME StokesGMRES
        FILES  StokesGMRES.cpp
        DEPENDS hyteg)
