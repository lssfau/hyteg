target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    FEFunctionWriter.hpp
    KeyValueStore.hpp
    SQL.cpp
    SQL.hpp
    Table.hpp
    TimingOutput.hpp
)

if( HYTEG_BUILD_WITH_ADIOS2 )
  add_subdirectory(ADIOS2)
endif()

add_subdirectory( VTKOutput )
