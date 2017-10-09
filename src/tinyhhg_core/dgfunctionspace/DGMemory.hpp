#pragma once

namespace hhg
{

inline uint_t DGVertexFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( level );
  return numDependencies * 2;
}

inline uint_t DGEdgeFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  const size_t num_cell_dofs = numDependencies * ( 2 * levelinfo::num_microedges_per_edge( level ) - 1 );
  return num_cell_dofs;
}

inline uint_t DGFaceFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( numDependencies );
  return levelinfo::num_microfaces_per_face( level );
}


}