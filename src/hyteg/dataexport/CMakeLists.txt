target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    KeyValueStore.hpp
    SQL.cpp
    SQL.hpp
    Table.hpp
    TimingOutput.hpp
)

add_subdirectory( VTKOutput )
