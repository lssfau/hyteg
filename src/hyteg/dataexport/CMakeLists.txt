target_sources( hyteg
    PRIVATE
    FEFunctionWriter.hpp
    SQL.cpp
    SQL.hpp
    TimingOutput.hpp
)

if( HYTEG_BUILD_WITH_ADIOS2 )
  add_subdirectory( ADIOS2 )
endif()

add_subdirectory( Gmsh )
add_subdirectory( LaTeX )
add_subdirectory( VTKOutput )
