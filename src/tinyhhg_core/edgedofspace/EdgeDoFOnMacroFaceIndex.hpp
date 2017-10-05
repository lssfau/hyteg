
#pragma once

#include "core/debug/Debug.h"

namespace hhg {
namespace EdgeDoFOnMacroFace {

using walberla::uint_t;

template<size_t Level>
constexpr inline uint_t indexFromVertex( const uint_t col, const uint_t row, const stencilDirection dir )
{
  typedef stencilDirection sD;

  const uint_t numMicroVerticesPerEdge = levelinfo::num_microvertices_per_edge( Level );
  const uint_t numMicroEdgesPerEdge    = levelinfo::num_microedges_per_edge( Level );

  const uint_t maxRowVertexIndex = numMicroVerticesPerEdge - col;
  const uint_t maxColVertexIndex = numMicroVerticesPerEdge - row;

  WALBERLA_ASSERT( col + row < numMicroVerticesPerEdge );
  WALBERLA_ASSERT( row <= maxRowVertexIndex );
  WALBERLA_ASSERT( col <= maxColVertexIndex );

  switch ( dir )
  {
  case sD::EDGE_HO_E:
    WALBERLA_ASSERT( row < maxRowVertexIndex );
    WALBERLA_ASSERT( col < maxColVertexIndex );
    return col + row * numMicroVerticesPerEdge - ( ( row + 1 ) * row ) / 2;
    break;
  case sD::EDGE_HO_W:
    WALBERLA_ASSERT( row <  maxRowVertexIndex );
    WALBERLA_ASSERT( col <= maxColVertexIndex );
    WALBERLA_ASSERT( col >  0 );
    return col - 1 + row * numMicroVerticesPerEdge - ( ( row + 1 ) * row ) / 2;
    break;
  default:
    WALBERLA_ASSERT( false );
    break;
  }
  return std::numeric_limits< uint_t >::max();
}

}
}
