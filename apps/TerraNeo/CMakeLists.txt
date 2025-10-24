add_subdirectory( Origin )
add_subdirectory( PlateVelocities ) 
add_subdirectory( SphericalHarmonics )

if( HYTEG_BUILD_WITH_ADIOS2 )
    add_subdirectory( TestApps )
    add_subdirectory( Utils )
endif()
