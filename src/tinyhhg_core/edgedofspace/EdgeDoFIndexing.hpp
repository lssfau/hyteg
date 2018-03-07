
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/levelinfo.hpp"

#include <cassert>

namespace hhg {
namespace edgedof {

constexpr uint_t levelToWidthAnyEdgeDoF( const uint_t & level )
{
  return levelinfo::num_microedges_per_edge( level );
}

constexpr uint_t levelToFaceSizeAnyEdgeDoF( const uint_t & level )
{
  return levelinfo::num_microedges_per_face( level ) / 3;
}

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

typedef stencilDirection sD;

/// Index of a horizontal edge DoF on a macro edge (only access to owned DoFs, no ghost layers).
inline constexpr uint_t horizontalIndex( const uint_t & level, const uint_t & col )
{
  return ::hhg::indexing::macroEdgeIndex( levelToWidthAnyEdgeDoF( level ), col );
};

/// Index of a horizontal edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
inline constexpr uint_t horizontalIndex( const uint_t & level, const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) - 1 );
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );

  const uint_t offset = numHorizontalDoFsOnEdge + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + horizontalIndex( level, col );
};

/// Index of a vertical edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
inline constexpr uint_t verticalIndex( const uint_t & level, const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) - 1  );
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );

  const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer + numOtherTypeDoFsOnGhostLayer + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + col;
};

/// Index of a diagonal edge DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
inline constexpr uint_t diagonalIndex( const uint_t & level, const uint_t & col, const uint_t & neighbor )
{
  const uint_t numHorizontalDoFsOnEdge       = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );
  const uint_t numHorizontalDoFsOnGhostLayer = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) - 1 );
  const uint_t numOtherTypeDoFsOnGhostLayer  = ::hhg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );

  const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer + neighbor * (numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer);

  return offset + col;
};

// Stencil access functions

inline constexpr uint_t indexFromHorizontalEdge( const uint_t & level, const uint_t & col, const stencilDirection & dir )
{
  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_C:
    return horizontalIndex( level, col );
  case sD::EDGE_DI_N:
    return diagonalIndex( level, col, 1 );
  case sD::EDGE_DI_S:
    return diagonalIndex( level, col, 0 );
  case sD::EDGE_VE_NW:
      return verticalIndex( level, col, 1 );
  case sD::EDGE_VE_SE:
      return verticalIndex( level, col, 0 );
  default:
    // assert( false );
    return std::numeric_limits< uint_t >::max();
  }
}


inline constexpr uint_t indexFromVertex( const uint_t & level, const uint_t & col, const stencilDirection & dir )
{
  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_W:
    return horizontalIndex( level, col - 1 );
  case sD::EDGE_HO_E:
    return horizontalIndex( level, col );
  case sD::EDGE_HO_NW:
    return horizontalIndex( level, col - 1, 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex( level, col - 1, 0 );
  case sD::EDGE_DI_SW:
    return diagonalIndex( level, col - 1, 0 );
  case sD::EDGE_DI_SE:
    return diagonalIndex( level, col, 0 );
  case sD::EDGE_DI_NW:
    return diagonalIndex( level, col - 1, 1 );
  case sD::EDGE_DI_NE:
    return diagonalIndex( level, col, 1 );
  case sD::EDGE_VE_N:
    return verticalIndex( level, col, 1 );
  case sD::EDGE_VE_S:
    return verticalIndex( level, col - 1, 0 );
  case sD::EDGE_VE_NW:
    return verticalIndex( level, col - 1, 1 );
  case sD::EDGE_VE_SE:
    return verticalIndex( level, col, 0 );
  default:
    // assert( false );
    return std::numeric_limits< uint_t >::max();
  }
}


constexpr std::array<stencilDirection ,2> neighborsOnEdgeFromVertex = {{ sD::EDGE_HO_E, sD::EDGE_HO_W}};
constexpr std::array<stencilDirection ,5> neighborsOnSouthFaceFromVertex = {{ sD::EDGE_DI_SW, sD::EDGE_VE_S, sD::EDGE_HO_SE, sD::EDGE_DI_SE, sD::EDGE_VE_SE}};
constexpr std::array<stencilDirection ,5> neighborsOnNorthFaceFromVertex = {{ sD::EDGE_DI_NE, sD::EDGE_VE_N, sD::EDGE_HO_NW, sD::EDGE_DI_NW, sD::EDGE_VE_NW}};

constexpr std::array<stencilDirection ,1> neighborsOnEdgeFromHorizontalEdge = {{ sD::EDGE_HO_C }};
constexpr std::array<stencilDirection ,2> neighborsOnSouthFaceFromHorizontalEdge = {{ sD::EDGE_DI_S, sD::EDGE_VE_SE }};
constexpr std::array<stencilDirection ,2> neighborsOnNorthFaceFromHorizontalEdge = {{ sD::EDGE_DI_N, sD::EDGE_VE_NW }};

class Iterator : public indexing::EdgeIterator
{
public:
  Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
    EdgeIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter )
  {}
};


} // namespace macroedge

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

typedef stencilDirection sD;

inline constexpr uint_t horizontalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
  return indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), col, row );
};

inline constexpr uint_t verticalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
  return 2 * levelToFaceSizeAnyEdgeDoF( level ) +
  indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), col, row );
}

inline constexpr uint_t diagonalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
  return levelToFaceSizeAnyEdgeDoF( level ) +
  indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), col, row );
}

// Stencil access functions

inline constexpr uint_t indexFromHorizontalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  switch( dir )
  {
  case sD::EDGE_HO_C:
    return horizontalIndex( level, col, row );
  case sD::EDGE_DI_N:
    return diagonalIndex( level, col, row );
  case sD::EDGE_DI_S:
    return diagonalIndex( level, col, row - 1 );
  case sD::EDGE_VE_NW:
      return verticalIndex( level, col, row );
  case sD::EDGE_VE_SE:
      return verticalIndex( level, col + 1, row - 1 );
  default:
    // assert( false );
    return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,5> neighborsFromHorizontalEdge =
  {{ sD::EDGE_HO_C,
     sD::EDGE_DI_S, sD::EDGE_VE_SE,
     sD::EDGE_DI_N, sD::EDGE_VE_NW
   }};


inline constexpr uint_t indexFromDiagonalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  switch( dir )
  {
  case sD::EDGE_DI_C:
    return diagonalIndex( level, col, row );
  case sD::EDGE_HO_N:
    return horizontalIndex( level, col, row + 1 );
  case sD::EDGE_HO_S:
    return horizontalIndex( level, col, row );
  case sD::EDGE_VE_W:
      return verticalIndex( level, col, row );
  case sD::EDGE_VE_E:
      return verticalIndex( level, col + 1, row );
  default:
    // assert( false );
    return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,5> neighborsFromDiagonalEdge =
  {{ sD::EDGE_DI_C,
     sD::EDGE_HO_S, sD::EDGE_VE_E,
     sD::EDGE_HO_N, sD::EDGE_VE_W
   }};

inline constexpr uint_t indexFromVerticalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  switch( dir )
  {
  case sD::EDGE_VE_C:
    return verticalIndex( level, col, row );
  case sD::EDGE_HO_NW:
    return horizontalIndex( level, col - 1, row + 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex( level, col, row );
  case sD::EDGE_DI_W:
      return diagonalIndex( level, col - 1, row );
  case sD::EDGE_DI_E:
      return diagonalIndex( level, col, row );
  default:
    // assert( false );
    return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,5> neighborsFromVerticalEdge =
  {{ sD::EDGE_VE_C,
     sD::EDGE_HO_SE, sD::EDGE_DI_E,
     sD::EDGE_HO_NW, sD::EDGE_DI_W
   }};

inline constexpr uint_t indexFromVertex( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
{
  // first  neighbor == south
  // second neighbor == north

  switch( dir )
  {
  case sD::EDGE_HO_W:
    return horizontalIndex( level, col - 1, row );
  case sD::EDGE_HO_E:
    return horizontalIndex( level, col, row );
  case sD::EDGE_HO_NW:
    return horizontalIndex( level, col - 1, row + 1 );
  case sD::EDGE_HO_SE:
    return horizontalIndex( level, col, row - 1 );
  case sD::EDGE_DI_SW:
    return diagonalIndex( level, col - 1, row - 1 );
  case sD::EDGE_DI_SE:
    return diagonalIndex( level, col, row - 1 );
  case sD::EDGE_DI_NW:
    return diagonalIndex( level, col - 1, row );
  case sD::EDGE_DI_NE:
    return diagonalIndex( level, col, row );
  case sD::EDGE_VE_N:
    return verticalIndex( level, col, row );
  case sD::EDGE_VE_S:
    return verticalIndex( level, col, row - 1 );
  case sD::EDGE_VE_NW:
    return verticalIndex( level, col - 1, row );
  case sD::EDGE_VE_SE:
    return verticalIndex( level, col + 1, row - 1 );
  default:
    WALBERLA_ASSERT( false );
    return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,12> neighborsFromVertex =
  {{ sD::EDGE_VE_S, sD::EDGE_HO_SE, sD::EDGE_DI_SE, sD::EDGE_VE_SE,
     sD::EDGE_HO_E, sD::EDGE_DI_NE, sD::EDGE_VE_N, sD::EDGE_HO_NW,
     sD::EDGE_DI_NW, sD::EDGE_VE_NW, sD::EDGE_HO_W, sD::EDGE_DI_SW
   }};

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromHorizontalEdge(const stencilDirection dir){
  switch(dir) {
    case sD::EDGE_DI_S:
      return 0;
    case sD::EDGE_VE_SE:
      return 1;
    case sD::EDGE_DI_N:
      return 2;
    case sD::EDGE_VE_NW:
      return 3;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromDiagonalEdge(const stencilDirection dir){
  switch(dir) {
    case sD::EDGE_HO_S:
      return 0;
    case sD::EDGE_VE_E:
      return 1;
    case sD::EDGE_HO_N:
      return 2;
    case sD::EDGE_VE_W:
      return 3;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromVerticalEdge(const stencilDirection dir){
  switch(dir) {
    case sD::EDGE_HO_S:
      return 0;
    case sD::EDGE_VE_E:
      return 1;
    case sD::EDGE_HO_N:
      return 2;
    case sD::EDGE_VE_W:
      return 3;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

// Iterators

class Iterator : public indexing::FaceIterator
{
public:
  Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
    FaceIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter )
  {}
};

class BorderIterator : public indexing::FaceBorderIterator
{
public:
  BorderIterator( const uint_t & level, const indexing::FaceBorderDirection & direction, const uint_t & offsetToCenter = 0, const uint_t & offsetFromVertices = 0 ) :
    FaceBorderIterator( levelinfo::num_microedges_per_edge( level ), direction, offsetToCenter, offsetFromVertices )
  {}
};

} // namespace macroface


/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are chosen such that the edge dofs on the south face from a macro edge are the first seven entries
/// otherwise the returned index would be out of bounds in the stencil memory array
constexpr inline uint_t stencilIndexFromVertex(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::EDGE_HO_W:
      return 0;
    case sD::EDGE_DI_SW:
      return 1;
    case sD::EDGE_VE_S:
      return 2;
    case sD::EDGE_HO_SE:
      return 3;
    case sD::EDGE_DI_SE:
      return 4;
    case sD::EDGE_VE_SE:
      return 5;
    case sD::EDGE_HO_E:
      return 6;
    case sD::EDGE_DI_NE:
      return 7;
    case sD::EDGE_VE_N:
      return 8;
    case sD::EDGE_HO_NW:
      return 9;
    case sD::EDGE_DI_NW:
      return 10;
    case sD::EDGE_VE_NW:
      return 11;
    default:
      return std::numeric_limits<size_t>::max();
  }
}


constexpr inline uint_t stencilIndexFromHorizontalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::EDGE_HO_C:
      return 0;
    case sD::EDGE_DI_S:
      return 1;
    case sD::EDGE_VE_SE:
      return 2;
    case sD::EDGE_DI_N:
      return 3;
    case sD::EDGE_VE_NW:
      return 4;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromDiagonalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::EDGE_DI_C:
      return 5;
    case sD::EDGE_HO_S:
      return 6;
    case sD::EDGE_VE_E:
      return 7;
    case sD::EDGE_HO_N:
      return 8;
    case sD::EDGE_VE_W:
      return 9;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromVerticalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::EDGE_VE_C:
      return 10;
    case sD::EDGE_HO_SE:
      return 11;
    case sD::EDGE_DI_E:
      return 12;
    case sD::EDGE_HO_NW:
      return 13;
    case sD::EDGE_DI_W:
      return 14;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

} // namespace edgedof
} // namespace hhg
