
#pragma once

#include "core/debug/Debug.h"

namespace hhg {
namespace EdgeDoFOnMacroEdge {

using walberla::uint_t;

template< size_t Level >
constexpr inline uint_t indexFromVertex( const uint_t pos, const stencilDirection dir )
{
  typedef stencilDirection sD;

  const uint_t numMicroVerticesPerEdge = levelinfo::num_microvertices_per_edge( Level );
  const uint_t numMicroEdgesPerEdge    = levelinfo::num_microedges_per_edge( Level );

  const uint_t maxVertexPos = numMicroVerticesPerEdge - 1;

  const uint_t numHorizontalEdgeDoFsCenter = numMicroEdgesPerEdge;
  const uint_t numHorizontalEdgeDoFsHalo   = numMicroEdgesPerEdge - 1;

  const uint_t numVerticalEdgeDoFsHalo     = numMicroEdgesPerEdge;
  const uint_t numDiagonalEdgeDoFsHalo     = numMicroEdgesPerEdge;

  // Only allowing micro-vertex positions that lie directly on the edge
  WALBERLA_ASSERT( pos < numMicroVerticesPerEdge );

  switch ( dir )
  {
  case sD::EDGE_HO_E:
    WALBERLA_ASSERT( pos < maxVertexPos );
    return pos;
    break;
  case sD::EDGE_HO_W:
    WALBERLA_ASSERT( pos <= maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1;
    break;
  case sD::EDGE_HO_NW:
    WALBERLA_ASSERT( pos <  maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter + numHorizontalEdgeDoFsHalo + numDiagonalEdgeDoFsHalo + numVerticalEdgeDoFsHalo;
    break;
  case sD::EDGE_HO_SE:
    WALBERLA_ASSERT( pos <  maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter;
    break;
  case sD::EDGE_VE_N:
    WALBERLA_ASSERT( pos < maxVertexPos );
    return pos + numHorizontalEdgeDoFsCenter + 2 * numHorizontalEdgeDoFsHalo + 2 * numDiagonalEdgeDoFsHalo + numVerticalEdgeDoFsHalo;
    break;
  case sD::EDGE_VE_S:
    WALBERLA_ASSERT( pos <= maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter + numHorizontalEdgeDoFsHalo + numDiagonalEdgeDoFsHalo;
    break;
  case sD::EDGE_VE_NW:
    WALBERLA_ASSERT( pos <= maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter + 2 * numHorizontalEdgeDoFsHalo + 2 * numDiagonalEdgeDoFsHalo + numVerticalEdgeDoFsHalo;
    break;
  case sD::EDGE_VE_SE:
    WALBERLA_ASSERT( pos < maxVertexPos );
    return pos + numHorizontalEdgeDoFsCenter + numHorizontalEdgeDoFsHalo + numDiagonalEdgeDoFsHalo;
    break;
  case sD::EDGE_DI_NW:
    WALBERLA_ASSERT( pos <= maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter + 2 * numHorizontalEdgeDoFsHalo + numDiagonalEdgeDoFsHalo + numVerticalEdgeDoFsHalo;
    break;
  case sD::EDGE_DI_NE:
    WALBERLA_ASSERT( pos < maxVertexPos );
    return pos + numHorizontalEdgeDoFsCenter + 2 * numHorizontalEdgeDoFsHalo + numDiagonalEdgeDoFsHalo + numVerticalEdgeDoFsHalo;
    break;
  case sD::EDGE_DI_SW:
    WALBERLA_ASSERT( pos <= maxVertexPos );
    WALBERLA_ASSERT( pos >  0 );
    return pos - 1 + numHorizontalEdgeDoFsCenter + numHorizontalEdgeDoFsHalo;
    break;
  case sD::EDGE_DI_SE:
    WALBERLA_ASSERT( pos < maxVertexPos );
    return pos + numHorizontalEdgeDoFsCenter + numHorizontalEdgeDoFsHalo;
    break;
  default:
    WALBERLA_ASSERT( false );
    break;
  }
  return std::numeric_limits< uint_t >::max();
}


template< size_t Level >
constexpr inline uint_t indexFromHorizontalEdge( const uint_t pos, const stencilDirection dir )
{
  typedef stencilDirection sD;

  WALBERLA_ASSERT( pos < levelinfo::num_microedges_per_edge( Level ) );

  switch ( dir )
  {
  case sD::EDGE_DI_N:
    return indexFromVertex< Level >( pos, sD::EDGE_DI_NE );
    break;
  case sD::EDGE_DI_S:
    return indexFromVertex< Level >( pos, sD::EDGE_DI_SE );
    break;
  case sD::EDGE_VE_NW:
    return indexFromVertex< Level >( pos, sD::EDGE_VE_N );
    break;
  case sD::EDGE_VE_SE:
    return indexFromVertex< Level >( pos, sD::EDGE_VE_SE );
    break;
  default:
    WALBERLA_ASSERT( false );
    break;
  }
  return std::numeric_limits< uint_t >::max();
}

}
}
