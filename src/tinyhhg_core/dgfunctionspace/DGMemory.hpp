#pragma once

#include "tinyhhg_core/Levelinfo.hpp"
#include "core/debug/CheckFunctions.h"

namespace hyteg {

using walberla::uint_t;

/// Vertex Memory layout
/// the vertex memory has two entries for each adjacent face,
/// where the first entry is the Gray Face DoF and the Second one is the Blue Face DoF
/// the Gray Face DoF is owned by the Vertex and the Blue Face DoF is owned by the Face
inline uint_t DGVertexFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return primitive.getNumNeighborFaces() * 2;
}

inline uint_t DGEdgeFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  const size_t num_cell_dofs = primitive.getNumNeighborFaces() * ( 2 * levelinfo::num_microedges_per_edge( level ) - 1 );
  return num_cell_dofs;
}

inline uint_t DGFaceFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return levelinfo::num_microfaces_per_face( level );
}


}
