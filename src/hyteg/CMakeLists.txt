add_library(hyteg)
target_link_libraries(hyteg PUBLIC core vtk)
if (HYTEG_BUILD_WITH_PETSC)
   target_link_libraries(hyteg PUBLIC PETSc::PETSc)
endif ()

if ( HYTEG_BUILD_WITH_ADIOS2 )
  if ( WALBERLA_BUILD_WITH_MPI )
    target_link_libraries( hyteg PUBLIC adios2::cxx11_mpi )
  else()
    target_link_libraries( hyteg PUBLIC adios2::cxx11 )
  endif()
endif()

set_target_properties(hyteg PROPERTIES CXX_EXTENSIONS OFF)
# There seem to be spurious warnings in gcc 12.1: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100366
# TODO: check if this is still the case for gcc 12.3
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12.1.0)
   target_compile_options(hyteg PUBLIC "-Wno-nonnull")
endif ()
target_sources(hyteg
      PRIVATE
      Algorithms.hpp
      BuildInfo.in.hpp
      Git.in.hpp
      HytegDefinitions.in.hpp
      Levelinfo.hpp
      LikwidWrapper.hpp
      Macros.hpp
      Math.hpp
      MeshQuality.hpp
      OpenMPManager.hpp
      PrimitiveID.hpp
      Stencil.hpp
      StencilDirections.hpp
      )

add_subdirectory( adaptiverefinement )
add_subdirectory( boundary )
add_subdirectory( checkpointrestore )
add_subdirectory( communication )
add_subdirectory( composites )
add_subdirectory( dataexport )
add_subdirectory( dg1functionspace )
add_subdirectory( dgfunctionspace )
add_subdirectory( edgedofspace )
add_subdirectory( egfunctionspace )
add_subdirectory( eigen )
add_subdirectory( elementwiseoperators )
add_subdirectory( fenics )
add_subdirectory( forms )
add_subdirectory( functions )
add_subdirectory( geometry )
add_subdirectory( gridtransferoperators )
add_subdirectory( indexing )
add_subdirectory( memory )
add_subdirectory( mesh )
add_subdirectory( misc )
add_subdirectory( mixedoperators )
add_subdirectory( n1e1functionspace )
add_subdirectory( numerictools )
add_subdirectory( operators )
add_subdirectory( p0functionspace )
add_subdirectory( p1functionspace )
add_subdirectory( p2functionspace )
add_subdirectory( petsc )
add_subdirectory( polynomial )
add_subdirectory( primitivedata )
add_subdirectory( primitives )
add_subdirectory( primitivestorage )
add_subdirectory( solvers )
add_subdirectory( sparseassembly )
add_subdirectory( trilinos )
add_subdirectory( types )
add_subdirectory( volumedofspace )
