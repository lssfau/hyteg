add_library( hyteg )
target_link_libraries( hyteg PUBLIC core vtk )
target_compile_features(hyteg PUBLIC cxx_std_17)
set_target_properties(hyteg PROPERTIES CXX_EXTENSIONS OFF )
target_sources( hyteg
      PRIVATE
      HytegDefinitions.in.hpp
      Git.in.hpp
      Levelinfo.hpp
      Macros.hpp
      OpenMPManager.hpp
      MeshQuality.hpp
      Math.hpp
      Stencil.hpp
      LikwidWrapper.hpp
      StencilDirections.hpp
      Algorithms.hpp
      PrimitiveID.hpp
      BuildInfo.in.hpp
      )

add_subdirectory( p2functionspace )
add_subdirectory( numerictools )
add_subdirectory( mesh )
add_subdirectory( edgedofspace )
add_subdirectory( misc )
add_subdirectory( adaptiverefinement )
add_subdirectory( gridtransferoperators )
add_subdirectory( types )
add_subdirectory( forms )
add_subdirectory( sparseassembly )
add_subdirectory( memory )
add_subdirectory( operators )
add_subdirectory( boundary )
add_subdirectory( composites )
add_subdirectory( solvers )
add_subdirectory( dgfunctionspace )
add_subdirectory( mixedoperators )
add_subdirectory( eigen )
add_subdirectory( elementwiseoperators )
add_subdirectory( fenics )
add_subdirectory( primitives )
add_subdirectory( dataexport )
add_subdirectory( primitivestorage )
add_subdirectory( celldofspace )
add_subdirectory( primitivedata )
add_subdirectory( functions )
add_subdirectory( trilinos )
add_subdirectory( petsc )
add_subdirectory( geometry )
add_subdirectory( indexing )
add_subdirectory( communication )
add_subdirectory( p1functionspace )
add_subdirectory( facedofspace )
add_subdirectory( polynomial )
