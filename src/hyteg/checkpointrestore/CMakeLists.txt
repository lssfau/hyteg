target_sources( hyteg
    PRIVATE
    CheckpointExporter.hpp
    CheckpointImporter.hpp
)

if( HYTEG_BUILD_WITH_ADIOS2 )
  add_subdirectory( ADIOS2 )
endif()
