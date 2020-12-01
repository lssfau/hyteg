waLBerla_add_module( DEPENDS core vtk )
target_compile_features(hyteg PUBLIC cxx_std_17)
set_target_properties(hyteg PROPERTIES CXX_EXTENSIONS OFF)

