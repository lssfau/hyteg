target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
)

add_subdirectory( VTKOutput )
