
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
inline constexpr uint_t horizontalIndex( const uint_t & col )
{
  return ::hhg::indexing::macroEdgeIndex< levelToWidthAnyEdgeDoF< level > >( col );
};

/// Index of a horizontal edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
template< uint_t level >
inline constexpr uint_t horizontalIndex( const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > - 1 >();
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();

  const uint_t offset = numHorizontalDoFsOnEdge + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + horizontalIndex< levelToWidthAnyEdgeDoF< level > >( col );
};

/// Index of a vertical edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
template< uint_t level >
inline constexpr uint_t verticalIndex( const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > - 1 >();
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();

  const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer + numOtherTypeDoFsOnGhostLayer + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + col;
};

/// Index of a diagonal edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
template< uint_t level >
inline constexpr uint_t diagonalIndex( const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > - 1 >();
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize< levelToWidthAnyEdgeDoF< level > >();

  const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + col;
};

// Stencil access functions

template< uint_t level >
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & col, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return horizontalIndex< level >( col );
  case sD::EDGE_DI_N:
    return diagonalIndex< level >( col, 1 );
  case sD::EDGE_DI_S:
    return diagonalIndex< level >( col, 0 );
  case sD::EDGE_VE_NW:
      return verticalIndex< level >( col, 1 );
  case sD::EDGE_VE_SE:
      return verticalIndex< level >( col, 0 );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}


template< uint_t level >
inline constexpr uint_t indexFromVertex( const uint_t & col, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_W:
    return horizontalIndex< level >( col - 1 );
  case sD::EDGE_HO_E:
    return horizontalIndex< level >( col     );
  case sD::EDGE_HO_NW:
    return horizontalIndex< level >( col - 1, 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex< level >( col - 1, 0 );
  case sD::EDGE_DI_SW:
    return diagonalIndex< level >( col - 1, 0 );
  case sD::EDGE_DI_SE:
    return diagonalIndex< level >( col    , 0 );
  case sD::EDGE_DI_NW:
    return diagonalIndex< level >( col - 1, 1 );
  case sD::EDGE_DI_NE:
    return diagonalIndex< level >( col    , 1 );
  case sD::EDGE_VE_N:
    return verticalIndex< level >( col    , 1 );
  case sD::EDGE_VE_S:
    return verticalIndex< level >( col - 1, 0 );
  case sD::EDGE_VE_NW:
    return verticalIndex< level >( col - 1, 1 );
  case sD::EDGE_VE_SE:
    return verticalIndex< level >( col    , 0 );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

} // namespace macroedge

// ##################
// ### Macro Face ###
// ##################

namespace macroface {


/// Direct access functions

template< uint_t level >
inline constexpr uint_t horizontalIndex( const uint_t & col, const uint_t & row )
{
  return macroFaceIndex< levelToWidthAnyEdgeDoF< level > >( col, row );
};

template< uint_t level >
inline constexpr uint_t verticalIndex( const uint_t & col, const uint_t & row )
{
  return 2 * levelToFaceSizeAnyEdgeDoF< level > + macroFaceIndex< levelToWidthAnyEdgeDoF< level > >( col, row );
}

template< uint_t level >
inline constexpr uint_t diagonalIndex( const uint_t & col, const uint_t & row )
{
  return levelToFaceSizeAnyEdgeDoF< level > + macroFaceIndex< levelToWidthAnyEdgeDoF< level > >( col, row );
}

// Stencil access functions

template< uint_t level >
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return horizontalIndex< level >( col, row );
  case sD::EDGE_DI_N:
    return diagonalIndex< level >( col, row     );
  case sD::EDGE_DI_S:
    return diagonalIndex< level >( col, row - 1 );
  case sD::EDGE_VE_NW:
      return verticalIndex< level >( col    , row     );
  case sD::EDGE_VE_SE:
      return verticalIndex< level >( col + 1, row - 1 );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

template< uint_t level >
inline constexpr uint_t indexFromDiagonalEdge( const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_DI_C:
    return diagonalIndex< level >( col, row );
  case sD::EDGE_HO_N:
    return horizontalIndex< level >( col, row + 1 );
  case sD::EDGE_HO_S:
    return horizontalIndex< level >( col, row     );
  case sD::EDGE_VE_W:
      return verticalIndex< level >( col    , row );
  case sD::EDGE_VE_E:
      return verticalIndex< level >( col + 1, row );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

template< uint_t level >
inline constexpr uint_t indexFromVerticalEdge( const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::EDGE_VE_C:
    return verticalIndex< level >( col, row );
  case sD::EDGE_HO_NW:
    return horizontalIndex< level >( col - 1, row + 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex< level >( col    , row     );
  case sD::EDGE_DI_W:
      return diagonalIndex< level >( col - 1, row );
  case sD::EDGE_DI_E:
      return diagonalIndex< level >( col    , row );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

template< uint_t level >
inline constexpr uint_t indexFromVertex( const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_W:
    return horizontalIndex< level >( col - 1, row     );
  case sD::EDGE_HO_E:
    return horizontalIndex< level >( col    , row     );
  case sD::EDGE_HO_NW:
    return horizontalIndex< level >( col - 1, row + 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex< level >( col    , row - 1 );
  case sD::EDGE_DI_SW:
    return diagonalIndex< level >( col - 1, row - 1 );
  case sD::EDGE_DI_SE:
    return diagonalIndex< level >( col    , row - 1 );
  case sD::EDGE_DI_NW:
    return diagonalIndex< level >( col - 1, row     );
  case sD::EDGE_DI_NE:
    return diagonalIndex< level >( col    , row     );
  case sD::EDGE_VE_N:
    return verticalIndex< level >( col    , row     );
  case sD::EDGE_VE_S:
    return verticalIndex< level >( col    , row - 1 );
  case sD::EDGE_VE_NW:
    return verticalIndex< level >( col - 1, row     );
  case sD::EDGE_VE_SE:
    return verticalIndex< level >( col + 1, row - 1 );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

// Iterators

template< uint_t level, uint_t offsetToCenter = 0 >
using Iterator = FaceIterator< levelToWidthAnyEdgeDoF< level >, offsetToCenter >;

template< uint_t level, uint_t offsetToCenter = 0 >
using BorderIterator = FaceBorderIterator< levelToWidthAnyEdgeDoF< level >, offsetToCenter >;

} // namespace macroface

} // namespace edgedof
} // namespace indexing
} // namespace hhg
