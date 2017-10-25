
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg {
namespace indexing {
namespace edgedof {

template< uint_t level >
constexpr uint_t levelToWidthAnyEdgeDoF = levelinfo::num_microedges_per_edge( level );

template< uint_t level >
constexpr uint_t levelToFaceSizeAnyEdgeDoF = levelinfo::num_microedges_per_face( level ) / 3;

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

/// Index of a horizontal edge DoF on a macro edge (only access to owned DoFs, no ghost layers).
template< uint_t level >
constexpr uint_t horizontalIndex( const uint_t & col )
{
  return ::hhg::indexing::macroEdgeIndex< levelToWidthAnyEdgeDoF< level > >( col );
};

/// Index of a horizontal edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
template< uint_t level >
constexpr uint_t horizontalIndex( const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > - 1 >();
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();

  const uint_t offset = numHorizontalDoFsOnEdge + neighbor * (numHorizontalDoFsOnGhostLayer + numOtherTypeDoFsOnGhostLayer);

  return offset + horizontalIndex< levelToWidthAnyEdgeDoF< level > >( col );
};

}

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

template< uint_t level >
constexpr uint_t horizontalIndex( const uint_t & col, const uint_t & row )
{
  return macroFaceIndex< levelToWidthAnyEdgeDoF< level > >( col, row );
};

// Stencil access functions



template< uint_t level >
constexpr uint_t indexFromVertex( const uint_t & col, const uint_t & row,
                                             const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return horizontalIndex< level >( col, row );

    // ...

  default:
    return 0;
    break;
  }
}

}

}
}
}
