set( WALBERLA_DIR    WALBERLA_DIR-NOTFOUND   CACHE  PATH  "waLBerla path"  )

if ( WALBERLA_DIR )
    # WALBERLA_DIR has to point to the waLBerla source directory
    # this command builds waLBerla (again) in the current build directory in the subfolder "walberla" (second argument)
    add_subdirectory( ${WALBERLA_DIR} walberla  )
    
    waLBerla_import()
    # Adds the 'src' and 'tests' directory of current app
    list( APPEND WALBERLA_MODULE_DIRS "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/tests" )
    list( REMOVE_DUPLICATES  WALBERLA_MODULE_DIRS )
    set ( WALBERLA_MODULE_DIRS  ${WALBERLA_MODULE_DIRS} CACHE INTERNAL "All folders that contain modules or tests" )
else()
    message( FATAL_ERROR "waLBerla not found - Use 'cmake -DWALBERLA_DIR=path_to_waLBerla_sources  pathToApplicationSources' "  )
endif()
    
