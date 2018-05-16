
#pragma once

#include "core/Abort.h"

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/levelinfo.hpp"

#include <cassert>

namespace hhg {
namespace vertexdof {

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

/// Index of a vertex DoF on a macro edge (only access to owned DoFs, no ghost layers).
inline constexpr uint_t index( const uint_t & level, const uint_t & x )
{
  return hhg::indexing::macroEdgeIndex( levelinfo::num_microvertices_per_edge( level ), x );
};

/// Index of a vertex DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor data, 1 to access second neighbor, ...
inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & neighbor )
{
  return hhg::indexing::macroEdgeSize( levelinfo::num_microvertices_per_edge( level ) )
         + neighbor * hhg::indexing::macroEdgeSize( levelinfo::num_microvertices_per_edge( level ) - 1 )
         + hhg::indexing::macroEdgeIndex( levelinfo::num_microvertices_per_edge( level ) - 1, x );
};

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
inline constexpr uint_t indexFromVertex( const uint_t & level, const uint_t & x, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_C:
      return index( level, x );
    case sD::VERTEX_E:
      return index( level, x + 1 );
    case sD::VERTEX_W:
      return index( level, x - 1 );
    case sD::VERTEX_N:
      return index( level, x, 1 );
    case sD::VERTEX_S:
      return index( level, x - 1, 0 );
    case sD::VERTEX_NW:
      return index( level, x - 1, 1 );
    case sD::VERTEX_SE:
      return index( level, x, 0 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

/// Have a look into the documentation to understand the calculations here
/// The west vertices have the same col index as the horizonal edge
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & level, const uint_t & x, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_W:
      return index( level, x );
    case sD::VERTEX_E:
      return index( level, x + 1 );
    case sD::VERTEX_SE:
      return index( level, x, 0 );
    case sD::VERTEX_NW:
      return index( level, x, 1 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

/// neighbor arrays from vertex dof
constexpr std::array<stencilDirection, 7> neighborsWithCenter               = {{ hhg::stencilDirection::VERTEX_C,
                                                                                 hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE,
                                                                                 hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_N,
                                                                                 hhg::stencilDirection::VERTEX_NW, hhg::stencilDirection::VERTEX_W
                                                                              }};
constexpr std::array<stencilDirection, 2> neighborsOnEdgeFromVertexDoF      = {{ hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_W }};
constexpr std::array<stencilDirection, 2> neighborsOnSouthFaceFromVertexDoF = {{ hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE }};
constexpr std::array<stencilDirection, 2> neighborsOnNorthFaceFromVertexDoF = {{ hhg::stencilDirection::VERTEX_N, hhg::stencilDirection::VERTEX_NW }};

/// neighbor arrays need to connect vertex dof and edge dof
constexpr std::array<stencilDirection, 2> neighborsOnEdgeFromHorizontalEdgeDoF      = {{ hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_W}};
constexpr std::array<stencilDirection, 1> neighborsOnSouthFaceFromHorizontalEdgeDoF = {{ hhg::stencilDirection::VERTEX_SE}};
constexpr std::array<stencilDirection, 1> neighborsOnNorthFaceFromHorizontalEdgeDoF = {{ hhg::stencilDirection::VERTEX_NW}};


/// Iterator over a vertex DoF macro edge.
/// See \ref EdgeIterator for more information.
class Iterator : public hhg::indexing::EdgeIterator
{
public:
  Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
    EdgeIterator( levelinfo::num_microvertices_per_edge( level ), offsetToCenter )
  {}
};

}

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

/// Index of a vertex DoF on a macro face (only access to owned DoFs, no ghost layers).
inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & y )
{
  return hhg::indexing::macroFaceIndex( levelinfo::num_microvertices_per_edge( level ), x, y );
};

/// Index of a vertex DoF on a ghost layer of a macro face.
/// \param neighbor 0 or 1 for the respective neighbor
inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & neighbor )
{
  assert( neighbor <= 1 );

  return hhg::indexing::macroFaceSize( levelinfo::num_microvertices_per_edge( level ) )
         + neighbor * hhg::indexing::macroFaceSize( levelinfo::num_microvertices_per_edge( level ) - 1 )
         + hhg::indexing::macroFaceIndex( levelinfo::num_microvertices_per_edge( level ) - 1, x, y );

};

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
inline constexpr uint_t indexFromVertex( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return index( level, x, y );
  case sD::VERTEX_E:
    return index( level, x + 1, y );
  case sD::VERTEX_W:
    return index( level, x - 1, y );
  case sD::VERTEX_N:
    return index( level, x, y + 1 );
  case sD::VERTEX_S:
    return index( level, x, y - 1 );
  case sD::VERTEX_NW:
    return index( level, x - 1, y + 1 );
  case sD::VERTEX_SE:
    return index( level, x + 1, y - 1 );
  case sD::VERTEX_TC:
    return index( level, x, y, 0 );
  case sD::VERTEX_TW:
    return index( level, x - 1, y, 0 );
  case sD::VERTEX_TS:
    return index( level, x, y - 1, 0 );
  case sD::VERTEX_TSE:
    return index( level, x - 1, y - 1, 0 );
  case sD::VERTEX_BC:
    return index( level, x - 1, y - 1, 1 );
  case sD::VERTEX_BN:
    return index( level, x - 1, y, 1 );
  case sD::VERTEX_BE:
    return index( level, x, y - 1, 1 );
  case sD::VERTEX_BNW:
    return index( level, x, y, 1 );
  default:
    return std::numeric_limits< uint_t >::max();
  }
}

/// Have a look into the documentation to understand the calculations here
/// The west vertex has the same col and row index as the horizonal edge
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_W:
      return index( level, x, y );
    case sD::VERTEX_E:
      return index( level, x + 1, y );
    case sD::VERTEX_SE:
      return index( level, x + 1, y - 1 );
    case sD::VERTEX_NW:
      return index( level, x, y + 1 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,4> neighborsFromHorizontalEdge =
  {{ stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E,
     stencilDirection::VERTEX_NW, stencilDirection::VERTEX_W
   }};

/// Have a look into the documentation to understand the calculations here
/// The south west vertex has the same col and row index as the horizonal edge
inline constexpr uint_t indexFromDiagonalEdge( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_SE:
      return index( level, x + 1, y );
    case sD::VERTEX_NE:
      return index( level, x + 1, y + 1 );
    case sD::VERTEX_NW:
      return index( level, x, y + 1 );
    case sD::VERTEX_SW:
      return index( level, x, y );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection ,4> neighborsFromDiagonalEdge =
  {{ stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NE,
     stencilDirection::VERTEX_NW, stencilDirection::VERTEX_SW
   }};

/// Have a look into the documentation to understand the calculations here
/// The south vertex has the same col and row index as the horizonal edge
inline constexpr uint_t indexFromVerticalEdge( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_S:
      return index( level, x, y );
    case sD::VERTEX_SE:
      return index( level, x + 1, y );
    case sD::VERTEX_N:
      return index( level, x, y + 1 );
    case sD::VERTEX_NW:
      return index( level, x - 1, y + 1 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

inline constexpr uint_t indexFromGrayFace( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_SW:
    return index( level, x, y );
  case sD::VERTEX_SE:
    return index( level, x + 1, y );
  case sD::VERTEX_NW:
    return index( level, x, y + 1 );
  default:
    return std::numeric_limits< uint_t >::max();
  }
}

inline constexpr uint_t indexFromBlueFace( const uint_t & level, const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_SE:
    return index( level, x + 1, y );
  case sD::VERTEX_NW:
    return index( level, x, y + 1 );
  case sD::VERTEX_NE:
    return index( level, x + 1, y + 1 );
  default:
    return std::numeric_limits< uint_t >::max();
  }
}

constexpr std::array<stencilDirection, 7> neighborsWithCenter     = {{ hhg::stencilDirection::VERTEX_C,
                                                                       hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE,
                                                                       hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_N,
                                                                       hhg::stencilDirection::VERTEX_NW, hhg::stencilDirection::VERTEX_W
                                                                     }};
constexpr std::array< stencilDirection, 6 > neighborsWithoutCenter = {{ hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE,
                                                                        hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_N,
                                                                        hhg::stencilDirection::VERTEX_NW, hhg::stencilDirection::VERTEX_W
                                                                      }};

constexpr std::array<stencilDirection ,4> neighborsFromVerticalEdge =
  {{ stencilDirection::VERTEX_S, stencilDirection::VERTEX_SE,
     stencilDirection::VERTEX_N, stencilDirection::VERTEX_NW
   }};

constexpr std::array< stencilDirection, 3 > neighborsFromGrayFace = {{ stencilDirection::VERTEX_SW, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NW }};
constexpr std::array< stencilDirection, 3 > neighborsFromBlueFace = {{ stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_NE }};

// Iterators

/// Iterator over a vertex DoF macro face.
/// See \ref FaceIterator for more information.
class Iterator : public hhg::indexing::FaceIterator
{
public:
  Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
    FaceIterator( levelinfo::num_microvertices_per_edge( level ), offsetToCenter )
  {}
};

/// Iterator over the border of a vertex DoF macro face.
/// See \ref FaceBorderIterator for more information.
class BorderIterator : public hhg::indexing::FaceBorderIterator
{
public:
  BorderIterator( const uint_t & level, const hhg::indexing::FaceBorderDirection & direction, const uint_t & offsetToCenter = 0, const uint_t & offsetFromVertices = 0) :
    FaceBorderIterator( levelinfo::num_microvertices_per_edge( level ), direction, offsetToCenter, offsetFromVertices )
  {}
};

} /// namespace macroface


// ##################
// ### Macro Cell ###
// ##################

namespace macrocell {

/// Index of a vertex DoF on a macro cell.
inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return hhg::indexing::macroCellIndex( levelinfo::num_microvertices_per_edge( level ), x, y, z );
}

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
inline constexpr uint_t indexFromVertex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_C:
      return index( level, x, y, z );
    case sD::VERTEX_W:
      return index( level, x - 1, y, z );
    case sD::VERTEX_E:
      return index( level, x + 1, y, z );
    case sD::VERTEX_N:
      return index( level, x, y + 1, z );
    case sD::VERTEX_S:
      return index( level, x, y - 1, z );
    case sD::VERTEX_NW:
      return index( level, x - 1, y + 1, z );
    case sD::VERTEX_SE:
      return index( level, x + 1, y - 1, z );
    case sD::VERTEX_TC:
      return index( level, x, y, z + 1 );
    case sD::VERTEX_TW:
      return index( level, x - 1, y, z + 1 );
    case sD::VERTEX_TS:
      return index( level, x, y - 1, z + 1 );
    case sD::VERTEX_TSE:
      return index( level, x + 1, y - 1, z + 1 );
    case sD::VERTEX_BC:
      return index( level, x, y, z - 1 );
    case sD::VERTEX_BN:
      return index( level, x, y + 1, z - 1 );
    case sD::VERTEX_BE:
      return index( level, x + 1, y, z - 1 );
    case sD::VERTEX_BNW:
      return index( level, x - 1, y + 1, z - 1 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}


/// neighbor arrays from vertex dof
constexpr std::array<stencilDirection, 15> neighborsWithCenter = {{ hhg::stencilDirection::VERTEX_C,
                                                                    hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE,
                                                                    hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_N,
                                                                    hhg::stencilDirection::VERTEX_NW, hhg::stencilDirection::VERTEX_W,
                                                                    hhg::stencilDirection::VERTEX_TC, hhg::stencilDirection::VERTEX_TW,
                                                                    hhg::stencilDirection::VERTEX_TS, hhg::stencilDirection::VERTEX_TSE,
                                                                    hhg::stencilDirection::VERTEX_BC, hhg::stencilDirection::VERTEX_BN,
                                                                    hhg::stencilDirection::VERTEX_BE, hhg::stencilDirection::VERTEX_BNW,
                                                                 }};

constexpr std::array<stencilDirection, 14> neighborsWithoutCenter = {{
                                                                       hhg::stencilDirection::VERTEX_S, hhg::stencilDirection::VERTEX_SE,
                                                                       hhg::stencilDirection::VERTEX_E, hhg::stencilDirection::VERTEX_N,
                                                                       hhg::stencilDirection::VERTEX_NW, hhg::stencilDirection::VERTEX_W,
                                                                       hhg::stencilDirection::VERTEX_TC, hhg::stencilDirection::VERTEX_TW,
                                                                       hhg::stencilDirection::VERTEX_TS, hhg::stencilDirection::VERTEX_TSE,
                                                                       hhg::stencilDirection::VERTEX_BC, hhg::stencilDirection::VERTEX_BN,
                                                                       hhg::stencilDirection::VERTEX_BE, hhg::stencilDirection::VERTEX_BNW,
                                                                    }};

// Iterators

/// Iterator over a vertex DoF macro cell.
/// See \ref CellIterator for more information.
class Iterator : public hhg::indexing::CellIterator
{
public:
  Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
    CellIterator( levelinfo::num_microvertices_per_edge( level ), offsetToCenter )
  {}
};

/// Iterator over the borders (faces) of a macro-cell.
/// See \ref CellBorderIterator for more information.
class BorderIterator : public hhg::indexing::CellBorderIterator
{
public:
  BorderIterator( const uint_t & level, const uint_t & vertex0, const uint_t & vertex1,
                  const uint_t & vertex2, const uint_t & offsetToCenter = 0 ) :
    CellBorderIterator( levelinfo::num_microvertices_per_edge( level ), vertex0, vertex1, vertex2, offsetToCenter )
  {}
};

} // namespace macrocell


// ################
// ### Stencils ###
// ################

/// Returns the logical index offset from a micro-vertex resulting from moving in the passed stencil direction
inline indexing::IndexIncrement logicalIndexOffsetFromVertex( const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_C:
      return indexing::IndexIncrement( 0, 0, 0 );
    case sD::VERTEX_W:
      return indexing::IndexIncrement( -1, 0, 0 );
    case sD::VERTEX_E:
      return indexing::IndexIncrement( 1, 0, 0 );
    case sD::VERTEX_N:
      return indexing::IndexIncrement( 0, 1, 0 );
    case sD::VERTEX_S:
      return indexing::IndexIncrement( 0, -1, 0 );
    case sD::VERTEX_NW:
      return indexing::IndexIncrement( -1, 1, 0 );
    case sD::VERTEX_SE:
      return indexing::IndexIncrement( 1, -1, 0 );
    case sD::VERTEX_TC:
      return indexing::IndexIncrement( 0, 0, 1 );
    case sD::VERTEX_TW:
      return indexing::IndexIncrement( -1, 0, 1 );
    case sD::VERTEX_TS:
      return indexing::IndexIncrement( 0, -1, 1 );
    case sD::VERTEX_TSE:
      return indexing::IndexIncrement( 1, -1, 1 );
    case sD::VERTEX_BC:
      return indexing::IndexIncrement( 0, 0, -1 );
    case sD::VERTEX_BN:
      return indexing::IndexIncrement( 0, 1, -1 );
    case sD::VERTEX_BE:
      return indexing::IndexIncrement( 1, 0, -1 );
    case sD::VERTEX_BNW:
      return indexing::IndexIncrement( -1, 1, -1 );
    default:
      WALBERLA_ASSERT( false, "Invalid stencil direction" );
      return indexing::IndexIncrement( std::numeric_limits< int >::max(),
                                       std::numeric_limits< int >::max(),
                                       std::numeric_limits< int >::max() );
  }
}

constexpr inline uint_t stencilIndexFromVertex( const stencilDirection dir )
{
  typedef stencilDirection sD;
  switch(dir) {
    case sD::VERTEX_S:
      return 0;
    case sD::VERTEX_SE:
      return 1;
    case sD::VERTEX_W:
      return 2;
    case sD::VERTEX_C:
      return 3;
    case sD::VERTEX_E:
      return 4;
    case sD::VERTEX_NW:
      return 5;
    case sD::VERTEX_N:
      return 6;
    case sD::VERTEX_TC:
      return 7;
    case sD::VERTEX_TW:
      return 8;
    case sD::VERTEX_TS:
      return 9;
    case sD::VERTEX_TSE:
      return 10;
    case sD::VERTEX_BC:
      return 11;
    case sD::VERTEX_BN:
      return 12;
    case sD::VERTEX_BE:
      return 13;
    case sD::VERTEX_BNW:
      return 14;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromHorizontalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::VERTEX_E:
      return 0;
    case sD::VERTEX_W:
      return 1;
    case sD::VERTEX_SE:
      return 2;
    case sD::VERTEX_NW:
      return 3;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromDiagonalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::VERTEX_SE:
      return 4;
    case sD::VERTEX_NE:
      return 5;
    case sD::VERTEX_NW:
      return 6;
    case sD::VERTEX_SW:
      return 7;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromVerticalEdge(const stencilDirection dir){
  typedef stencilDirection sD;
  switch(dir) {
    case sD::VERTEX_S:
      return 8;
    case sD::VERTEX_SE:
      return 9;
    case sD::VERTEX_N:
      return 10;
    case sD::VERTEX_NW:
      return 11;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromGrayFace( const stencilDirection & dir )
{
  typedef stencilDirection sD;
  switch (dir) {
    case sD::VERTEX_SW:
      return 0;
    case sD::VERTEX_SE:
      return 1;
    case sD::VERTEX_NW:
      return 2;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr inline uint_t stencilIndexFromBlueFace( const stencilDirection & dir )
{
  typedef stencilDirection sD;
  switch (dir) {
    case sD::VERTEX_SE:
      return 0;
    case sD::VERTEX_NW:
      return 1;
    case sD::VERTEX_NE:
      return 2;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

inline bool isVertexOnBoundary(const uint_t level, const hhg::indexing::Index& idx){
  if (idx.row() == 0 ){
    return true;
  } else if( idx.col() == 0 ){
    return true;
  } else if( (idx.row() + idx.col()) == ( hhg::levelinfo::num_microvertices_per_edge( level ) - 1 ) ){
    return true;
  } else {
    return false;
  }
}

} /// namespace vertexdof
} /// namespace hhg
