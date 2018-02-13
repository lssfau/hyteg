
#pragma once

#include "core/Abort.h"

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/levelinfo.hpp"

namespace hhg {
namespace vertexdof {

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

/// Index of a vertex DoF on a macro edge (only access to owned DoFs, no ghost layers).
template< uint_t level >
inline constexpr uint_t index( const uint_t & x )
{
  return hhg::indexing::macroEdgeIndex< levelinfo::num_microvertices_per_edge( level ) >( x );
};

/// Index of a vertex DoF on a ghost layer of a macro edge.
/// \param neighbor 0 to access the first neighbor data, 1 to access second neighbor, ...
template< uint_t level >
inline constexpr uint_t index( const uint_t & x, const uint_t & neighbor )
{
  return              hhg::indexing::macroEdgeSize< levelinfo::num_microvertices_per_edge( level )     >()
         + neighbor * hhg::indexing::macroEdgeSize< levelinfo::num_microvertices_per_edge( level ) - 1 >()
         + hhg::indexing::macroEdgeIndex< levelinfo::num_microvertices_per_edge( level ) - 1 >( x );
};

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
template< uint_t level >
inline constexpr uint_t indexFromVertex( const uint_t & x, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_C:
      return index< level >( x    );
    case sD::VERTEX_E:
      return index< level >( x + 1);
    case sD::VERTEX_W:
      return index< level >( x - 1);
    case sD::VERTEX_N:
      return index< level >( x    , 1);
    case sD::VERTEX_S:
      return index< level >( x - 1, 0);
    case sD::VERTEX_NW:
      return index< level >( x - 1, 1);
    case sD::VERTEX_SE:
      return index< level >( x    , 0);
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

/// Have a look into the documentation to understand the calculations here
/// The west vertices have the same col index as the horizonal edge
template< uint_t level >
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & x, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_W:
      return index< level >( x    );
    case sD::VERTEX_E:
      return index< level >( x + 1);
    case sD::VERTEX_SE:
      return index< level >( x , 0);
    case sD::VERTEX_NW:
      return index< level >( x , 1);
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
template< uint_t level >
inline constexpr uint_t index( const uint_t & x, const uint_t & y )
{
  return hhg::indexing::macroFaceIndex< levelinfo::num_microvertices_per_edge( level ) >( x, y );
};

/// Index of a vertex DoF on a ghost layer of a macro face.
/// \param neighbor 0 or 1 for the respective neighbor
template< uint_t level >
inline constexpr uint_t index( const uint_t & x, const uint_t & y, const uint_t & neighbor )
{
  WALBERLA_ASSERT( neighbor <= 1 );

  return              hhg::indexing::macroFaceSize< levelinfo::num_microvertices_per_edge( level )     >()
         + neighbor * hhg::indexing::macroFaceSize< levelinfo::num_microvertices_per_edge( level ) - 1 >()
         + hhg::indexing::macroFaceIndex< levelinfo::num_microvertices_per_edge( level ) - 1 >( x, y );

};

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
template< uint_t level >
inline constexpr uint_t indexFromVertex( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_C:
    return index< level >( x    , y     );
  case sD::VERTEX_E:
    return index< level >( x + 1, y     );
  case sD::VERTEX_W:
    return index< level >( x - 1, y     );
  case sD::VERTEX_N:
    return index< level >( x    , y + 1 );
  case sD::VERTEX_S:
    return index< level >( x    , y - 1 );
  case sD::VERTEX_NW:
    return index< level >( x - 1, y + 1 );
  case sD::VERTEX_SE:
    return index< level >( x + 1, y - 1 );
  default:
    return std::numeric_limits< uint_t >::max();
  }
}

/// Have a look into the documentation to understand the calculations here
/// The west vertex has the same col and row index as the horizonal edge
template< uint_t level >
inline constexpr uint_t indexFromHorizontalEdge( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_W:
      return index< level >( x    , y);
    case sD::VERTEX_E:
      return index< level >( x + 1, y);
    case sD::VERTEX_SE:
      return index< level >( x + 1, y - 1);
    case sD::VERTEX_NW:
      return index< level >( x    , y + 1);
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
template< uint_t level >
inline constexpr uint_t indexFromDiagonalEdge( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_SE:
      return index< level >( x + 1, y);
    case sD::VERTEX_NE:
      return index< level >( x + 1, y + 1);
    case sD::VERTEX_NW:
      return index< level >( x    , y + 1);
    case sD::VERTEX_SW:
      return index< level >( x    , y    );
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
template< uint_t level >
inline constexpr uint_t indexFromVerticalEdge( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_S:
      return index< level >( x    , y    );
    case sD::VERTEX_SE:
      return index< level >( x + 1, y    );
    case sD::VERTEX_N:
      return index< level >( x    , y + 1);
    case sD::VERTEX_NW:
      return index< level >( x - 1, y + 1);
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

template< uint_t level >
inline constexpr uint_t indexFromGrayFace( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_SW:
    return index< level >( x    , y     );
  case sD::VERTEX_SE:
    return index< level >( x + 1, y     );
  case sD::VERTEX_NW:
    return index< level >( x    , y + 1 );
  default:
    return std::numeric_limits< uint_t >::max();
  }
}

template< uint_t level >
inline constexpr uint_t indexFromBlueFace( const uint_t & x, const uint_t & y, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
  case sD::VERTEX_SE:
    return index< level >( x + 1, y     );
  case sD::VERTEX_NW:
    return index< level >( x    , y + 1 );
  case sD::VERTEX_NE:
    return index< level >( x + 1, y + 1 );
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
  BorderIterator( const uint_t & level, const hhg::indexing::FaceBorderDirection & direction, const uint_t & offsetToCenter = 0 ) :
    FaceBorderIterator( levelinfo::num_microvertices_per_edge( level ), direction, offsetToCenter )
  {}
};

} /// namespace macroface


// ##################
// ### Macro Cell ###
// ##################

namespace macrocell {

/// Index of a vertex DoF on a macro cell.
template< uint_t level >
inline constexpr uint_t index( const uint_t & x, const uint_t & y, const uint_t & z )
{
  return hhg::indexing::macroCellIndex< levelinfo::num_microvertices_per_edge( level ) >( x, y, z );
}

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
template< uint_t level >
inline constexpr uint_t indexFromVertex( const uint_t & x, const uint_t & y, const uint_t & z, const stencilDirection & dir )
{
  typedef stencilDirection sD;

  switch( dir )
  {
    case sD::VERTEX_C:
      return index< level >( x    , y    , z     );
    case sD::VERTEX_W:
      return index< level >( x - 1, y    , z     );
    case sD::VERTEX_E:
      return index< level >( x + 1, y    , z     );
    case sD::VERTEX_N:
      return index< level >( x    , y + 1, z     );
    case sD::VERTEX_S:
      return index< level >( x    , y - 1, z     );
    case sD::VERTEX_NW:
      return index< level >( x - 1, y + 1, z     );
    case sD::VERTEX_SE:
      return index< level >( x + 1, y - 1, z     );
    case sD::VERTEX_BC:
      return index< level >( x    , y    , z + 1 );
    case sD::VERTEX_BW:
      return index< level >( x - 1, y    , z + 1 );
    case sD::VERTEX_BS:
      return index< level >( x    , y - 1, z + 1 );
    case sD::VERTEX_BSW:
      return index< level >( x - 1, y - 1, z + 1 );
    case sD::VERTEX_FC:
      return index< level >( x    , y    , z - 1 );
    case sD::VERTEX_FN:
      return index< level >( x    , y + 1, z - 1 );
    case sD::VERTEX_FE:
      return index< level >( x + 1, y    , z - 1 );
    case sD::VERTEX_FNE:
      return index< level >( x + 1, y + 1, z - 1 );
    default:
      return std::numeric_limits< uint_t >::max();
  }
}

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
  BorderIterator( const uint_t & level, const uint_t & vertex0, const uint_t & vertex1, const uint_t & vertex2 ) :
    CellBorderIterator( levelinfo::num_microvertices_per_edge( level ), vertex0, vertex1, vertex2 )
  {}
};

} // namespace macrocell


// ################
// ### Stencils ###
// ################

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

} /// namespace vertexdof
} /// namespace hhg
