waLBerla_add_executable( NAME 01_PrimitiveStorage
        FILES 01_PrimitiveStorage.cpp
        DEPENDS tinyhhg_core )

waLBerla_add_executable( NAME 02_PrimitiveData
        FILES 02_PrimitiveData.cpp
        DEPENDS tinyhhg_core )

waLBerla_add_executable( NAME 03_Communication
        FILES 03_Communication.cpp
        DEPENDS tinyhhg_core )

waLBerla_add_executable( NAME 04_Indexing
        FILES 04_Indexing.cpp
        DEPENDS tinyhhg_core )

waLBerla_add_executable( NAME FullAppP1GMG
        FILES 05_FullAppP1GMG/P1GMG.cpp
        DEPENDS tinyhhg_core core)
