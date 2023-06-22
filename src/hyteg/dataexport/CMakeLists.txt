target_sources( hyteg
    PRIVATE
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
)

add_subdirectory( VTKOutput )
