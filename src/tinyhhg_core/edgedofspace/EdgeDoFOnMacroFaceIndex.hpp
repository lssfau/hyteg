
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

  const uint_t numMicroVerticesInRow      = numMicroVerticesPerEdge - row;
  const uint_t numMicroEdgesInRow         = numMicroEdgesPerEdge - row;
  const uint_t numMicroEdgesInRowBelow    = numMicroEdgesPerEdge - (row - 1);

  const uint_t sumZeroToRow = ( ( row + 1 ) * row ) / 2;

  const uint_t numHorizontalEdgeDoFsBelowRow = row * numMicroVerticesPerEdge - sumZeroToRow;

  const uint_t numOverallHorizontalEdgeDoFs = (( numMicroEdgesPerEdge * numMicroEdgesPerEdge ) - numMicroEdgesPerEdge ) / 2 + numMicroEdgesPerEdge;
  const uint_t numOverallDiagonalEdgeDoFs   = numOverallHorizontalEdgeDoFs;
  const uint_t numOverallVerticalEdgeDoFs   = numOverallHorizontalEdgeDoFs;

  WALBERLA_ASSERT( col + row < numMicroVerticesPerEdge );

  switch ( dir )
  {
  case sD::EDGE_HO_E:
  {
    WALBERLA_ASSERT( col + row < numMicroVerticesPerEdge - 1 );
    const uint_t result = col + numHorizontalEdgeDoFsBelowRow;
    WALBERLA_ASSERT( result < numOverallHorizontalEdgeDoFs );
    return result;
    break;
  }
  case sD::EDGE_HO_W:
  {
    WALBERLA_ASSERT( col >  0 );
    const uint_t result = col - 1 + numHorizontalEdgeDoFsBelowRow;
    WALBERLA_ASSERT( result < numOverallHorizontalEdgeDoFs );
    return result;
    break;
  }
  case sD::EDGE_HO_NW:
  {
    WALBERLA_ASSERT( col > 0 );
    WALBERLA_ASSERT( col + row < numMicroVerticesPerEdge - 1 );
    const uint_t result = col - 1 + numHorizontalEdgeDoFsBelowRow + numMicroEdgesInRow;
    WALBERLA_ASSERT( result < numOverallHorizontalEdgeDoFs );
    return result;
    break;
  }
  case sD::EDGE_HO_SE:
  {
    WALBERLA_ASSERT( row > 0 );
    const uint_t result = col + numHorizontalEdgeDoFsBelowRow - numMicroEdgesInRowBelow;
    WALBERLA_ASSERT( result < numOverallHorizontalEdgeDoFs );
    return result;
    break;
  }
  case sD::EDGE_VE_N:
  {
    WALBERLA_ASSERT( col + row < numMicroVerticesPerEdge - 1 );
    const uint_t result = numOverallHorizontalEdgeDoFs + numOverallDiagonalEdgeDoFs + col + numHorizontalEdgeDoFsBelowRow;
    WALBERLA_ASSERT( result >= numOverallHorizontalEdgeDoFs + numOverallDiagonalEdgeDoFs );
    WALBERLA_ASSERT( result <  numOverallHorizontalEdgeDoFs + numOverallDiagonalEdgeDoFs + numOverallVerticalEdgeDoFs );
    return result;
    break;
  }
  default:
    WALBERLA_ASSERT( false );
    break;
  }
  return std::numeric_limits< uint_t >::max();
}

}
}
