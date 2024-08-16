waLBerla_add_executable( NAME BH.01_PrimitiveData
        FILES BH.01_PrimitiveData.cpp
        DEPENDS hyteg )
if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME BH.01_PrimitiveData )
endif()

waLBerla_add_executable( NAME BH.02_Communication
        FILES BH.02_Communication.cpp
        DEPENDS hyteg )
if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME BH.02_Communication )
endif()

waLBerla_add_executable( NAME BH.03_Indexing
        FILES BH.03_Indexing.cpp
        DEPENDS hyteg )
if(CMAKE_BUILD_TYPE MATCHES "Release")
   waLBerla_execute_test( NO_MODULE_LABEL NAME BH.03_Indexing )
endif()
