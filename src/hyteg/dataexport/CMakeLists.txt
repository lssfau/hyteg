target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    KeyValueStore.hpp
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
    Table.hpp
)

add_subdirectory( VTKOutput )
