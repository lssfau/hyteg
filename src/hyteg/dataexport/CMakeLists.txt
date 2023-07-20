target_sources( hyteg
    PRIVATE
    FEFunctionRegistry.hpp
    FEFunctionWriter.hpp
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
)

if( HYTEG_BUILD_WITH_ADIOS2 )
  add_subdirectory(ADIOS2)
endif()

add_subdirectory( LaTeX )
add_subdirectory( VTKOutput )
