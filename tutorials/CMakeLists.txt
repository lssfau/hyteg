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

add_subdirectory(05_FullAppP1GMG)
add_subdirectory(06_FullAppPlumeInCube)
add_subdirectory(07_IsoviscousConvectionAnnulus)
add_subdirectory(08_CahnHilliard)
add_subdirectory(09_BoundaryConditions)
add_subdirectory(10_DGAMR)
add_subdirectory(11_CouetteFlow)
add_subdirectory(12_UnresolvedParticles)
