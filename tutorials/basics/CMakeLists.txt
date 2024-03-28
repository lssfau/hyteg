waLBerla_add_executable( NAME 01_PrimitiveStorage
        FILES 01_PrimitiveStorage.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME 02_PrimitiveData
        FILES 02_PrimitiveData.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME 03_Communication
        FILES 03_Communication.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME 04_Indexing
        FILES 04_Indexing.cpp
        DEPENDS hyteg )

add_subdirectory(09_BoundaryConditions)
