waLBerla_link_files_to_builddir( data )

add_executable       ( adaptiveRefinement adaptiveRefinement.cpp )
target_link_libraries( adaptiveRefinement hyteg opgen-mass constant_stencil_operator )

add_executable       ( fullAMR fullAMR.cpp )
target_link_libraries( fullAMR hyteg opgen-mass constant_stencil_operator )