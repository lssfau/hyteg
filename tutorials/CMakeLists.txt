if(NOT EXISTS ${CMAKE_BINARY_DIR}/data/meshes)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/data/meshes")
endif()

waLBerla_link_files_to_builddir(../data/meshes/*.msh)

waLBerla_add_executable( NAME 01_PrimitiveStorage
        FILES 01_PrimitiveStorage.cpp
        DEPENDS tinyhhg_core )
        
waLBerla_add_executable( NAME 02_PrimitiveData
        FILES 02_PrimitiveData.cpp
        DEPENDS tinyhhg_core )