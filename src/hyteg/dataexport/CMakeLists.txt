target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    FEFunctionWriter.hpp
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
)

add_subdirectory( VTKOutput )
