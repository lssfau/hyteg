if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/output)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
endif()

waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable( NAME SPHdemo
        FILES SPHdemo.cpp
        DEPENDS hyteg)

if (NOT APPLE)
    waLBerla_add_executable(NAME EigFunLaplace
          FILES EigFunLaplace.cpp
          DEPENDS hyteg)
endif ()

waLBerla_add_executable( NAME bmark
        FILES bmark.cpp
        DEPENDS hyteg)

