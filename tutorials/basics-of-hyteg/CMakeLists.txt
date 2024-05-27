waLBerla_add_executable( NAME BH.01_PrimitiveData
        FILES BH.01_PrimitiveData.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME BH.02_Communication
        FILES BH.02_Communication.cpp
        DEPENDS hyteg )

waLBerla_add_executable( NAME BH.03_Indexing
        FILES BH.03_Indexing.cpp
        DEPENDS hyteg )
