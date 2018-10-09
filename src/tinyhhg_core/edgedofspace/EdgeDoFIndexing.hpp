
#pragma once

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"
#include "tinyhhg_core/indexing/MacroCellIndexing.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/Levelinfo.hpp"

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

enum class EdgeDoFOrientation : walberla::uint_t
{
  X,
  Y,
  Z,
  XY,
  XZ,
  YZ,
  XYZ,
  INVALID,
};

const std::array< EdgeDoFOrientation, 7 > allEdgeDoFOrientations = { EdgeDoFOrientation::X, EdgeDoFOrientation::Y, EdgeDoFOrientation::Z,
                                                                     EdgeDoFOrientation::XY, EdgeDoFOrientation::XZ, EdgeDoFOrientation::YZ,
                                                                     EdgeDoFOrientation::XYZ,};

const std::array< EdgeDoFOrientation, 3 > faceLocalEdgeDoFOrientations = { EdgeDoFOrientation::X, EdgeDoFOrientation::Y, EdgeDoFOrientation::XY };

const std::map< EdgeDoFOrientation, std::string > edgeDoFOrientationToString = {
  { EdgeDoFOrientation::X, "X"},
  { EdgeDoFOrientation::Y, "Y"},
  { EdgeDoFOrientation::Z, "Z"},
  { EdgeDoFOrientation::XY, "XY"},
  { EdgeDoFOrientation::XZ, "XZ"},
  { EdgeDoFOrientation::YZ, "YZ"},
  { EdgeDoFOrientation::XYZ, "XYZ"},
  { EdgeDoFOrientation::INVALID, "INVALID"},
};

inline std::ostream& operator<<(std::ostream& out, const EdgeDoFOrientation ornt)
{
   const std::string& orntString = edgeDoFOrientationToString.at(ornt);
   return out << orntString;
}

/// \brief Given two logical vertexdof indices, this function returns the appropriate edgedof orientation.
inline EdgeDoFOrientation calcEdgeDoFOrientation( const indexing::IndexIncrement & vertexIndex0, const indexing::IndexIncrement & vertexIndex1 )
{
  const indexing::IndexIncrement offset = vertexIndex1 - vertexIndex0;
  const uint_t x = std::abs( offset.x() );
  const uint_t y = std::abs( offset.y() );
  const uint_t z = std::abs( offset.z() );

  WALBERLA_ASSERT_GREATER( x + y + z, 0 );
  WALBERLA_ASSERT_LESS_EQUAL( x, 1 );
  WALBERLA_ASSERT_LESS_EQUAL( y, 1 );
  WALBERLA_ASSERT_LESS_EQUAL( z, 1 );

  if ( x == 1 && y == 0 && z == 0 )
    return EdgeDoFOrientation::X;
  if ( x == 0 && y == 1 && z == 0 )
    return EdgeDoFOrientation::Y;
  if ( x == 0 && y == 0 && z == 1 )
    return EdgeDoFOrientation::Z;
  if ( x == 1 && y == 1 && z == 0 )
    return EdgeDoFOrientation::XY;
  if ( x == 1 && y == 0 && z == 1 )
    return EdgeDoFOrientation::XZ;
  if ( x == 0 && y == 1 && z == 1 )
    return EdgeDoFOrientation::YZ;
  if ( x == 1 && y == 1 && z == 1 )
    return EdgeDoFOrientation::XYZ;

  WALBERLA_ASSERT( false, "Invalid index offset." );
  return EdgeDoFOrientation::INVALID;
}

/// \brief Given two logical vertexdof indices, this function returns the logical index of the edgedof inbetween.
/// This function also implicitly calculates the orientation.
inline indexing::IndexIncrement calcEdgeDoFIndex( const indexing::IndexIncrement & vertexIndex0, const indexing::IndexIncrement & vertexIndex1 )
{
  const EdgeDoFOrientation orientation = calcEdgeDoFOrientation( vertexIndex0, vertexIndex1 );
  switch ( orientation )
  {
  case EdgeDoFOrientation::X:
    return vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1;
  case EdgeDoFOrientation::Y:
    return vertexIndex0.y() < vertexIndex1.y() ? vertexIndex0 : vertexIndex1;
  case EdgeDoFOrientation::Z:
    return vertexIndex0.z() < vertexIndex1.z() ? vertexIndex0 : vertexIndex1;
  case EdgeDoFOrientation::XY:
    return (vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1) + indexing::IndexIncrement( 0, -1,  0 );
  case EdgeDoFOrientation::XZ:
    return (vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1) + indexing::IndexIncrement( 0,  0, -1 );
  case EdgeDoFOrientation::YZ:
    return (vertexIndex0.y() < vertexIndex1.y() ? vertexIndex0 : vertexIndex1) + indexing::IndexIncrement( 0,  0, -1 );
  case EdgeDoFOrientation::XYZ:
    return (vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1) + indexing::IndexIncrement( 0, -1,  0 );
  default:
    WALBERLA_ASSERT( false, "Invaild index offset." );
    return indexing::IndexIncrement( std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max() );
  }
}


/// \brief Given an orientation, this function returns the logical index offsets of the two neighboring vertexdofs on that edge.
/// If the resulting offsets are added to an logical edge index, the resulting logical indices are those of the neighboring vertices.
inline std::array< indexing::IndexIncrement, 2 > calcNeighboringVertexDoFIndices( const edgedof::EdgeDoFOrientation & orientation )
{
  std::array< indexing::IndexIncrement, 2 > vertexIndices = {
    indexing::IndexIncrement( std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max() ),
    indexing::IndexIncrement( std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max(), std::numeric_limits< uint_t >::max() )
  };
  vertexIndices[0] = indexing::IndexIncrement( 0, 0, 0 );

  switch ( orientation )
  {
    case EdgeDoFOrientation::X:
      vertexIndices[1] = indexing::IndexIncrement( 1, 0, 0 );
      break;
    case EdgeDoFOrientation::Y:
      vertexIndices[1] = indexing::IndexIncrement( 0, 1, 0 );
      break;
    case EdgeDoFOrientation::Z:
      vertexIndices[1] = indexing::IndexIncrement( 0, 0, 1 );
      break;
    case EdgeDoFOrientation::XY:
      vertexIndices[0] = indexing::IndexIncrement( 1, 0, 0 );
      vertexIndices[1] = indexing::IndexIncrement( 0, 1, 0 );
      break;
    case EdgeDoFOrientation::XZ:
      vertexIndices[0] = indexing::IndexIncrement( 1, 0, 0 );
      vertexIndices[1] = indexing::IndexIncrement( 0, 0, 1 );
      break;
    case EdgeDoFOrientation::YZ:
      vertexIndices[0] = indexing::IndexIncrement( 0, 1, 0 );
      vertexIndices[1] = indexing::IndexIncrement( 0, 0, 1 );
      break;
    case EdgeDoFOrientation::XYZ:
      vertexIndices[0] = indexing::IndexIncrement( 0, 1, 0 );
      vertexIndices[1] = indexing::IndexIncrement( 1, 0, 1 );
      break;
    default:
      WALBERLA_ASSERT( false, "Invaild index offset." );
  }

  return vertexIndices;
}

/// \brief Given to local vertex ids on the tetrahedron this function return the orientation of the edge inbetween those vertices
inline EdgeDoFOrientation getEdgeDoFOrientationFromLocalIDs( const uint_t & vertexIndex0, const uint_t & vertexIndex1 ){
   WALBERLA_ASSERT_UNEQUAL( vertexIndex0, vertexIndex1 );
   WALBERLA_ASSERT_LESS_EQUAL( vertexIndex0, 3 );
   WALBERLA_ASSERT_LESS_EQUAL( vertexIndex1, 3 );
   indexing::IndexIncrement vertexIndexII0, vertexIndexII1;
   switch ( vertexIndex0 ){
      case 0:
         vertexIndexII0 = indexing::IndexIncrement( 0, 0, 0 );
         break;
      case 1:
         vertexIndexII0 = indexing::IndexIncrement( 1, 0, 0 );
         break;
      case 2:
         vertexIndexII0 = indexing::IndexIncrement( 0, 1, 0 );
         break;
      case 3:
         vertexIndexII0 = indexing::IndexIncrement( 0, 0, 1 );
         break;
      default:
         WALBERLA_ABORT("Wrong vertex ID");

   }
   switch ( vertexIndex1 ){
      case 0:
         vertexIndexII1 = indexing::IndexIncrement( 0, 0, 0 );
         break;
      case 1:
         vertexIndexII1 = indexing::IndexIncrement( 1, 0, 0 );
         break;
      case 2:
         vertexIndexII1 = indexing::IndexIncrement( 0, 1, 0 );
         break;
      case 3:
         vertexIndexII1 = indexing::IndexIncrement( 0, 0, 1 );
         break;
      default:
         WALBERLA_ABORT("Wrong vertex ID");
   }

   return calcEdgeDoFOrientation(vertexIndexII0, vertexIndexII1);

}

/// \brief converts the edge orientation of the reference face into the correct orientation on the tetrahedron according to the position
/// and the rotation of the face.
inline EdgeDoFOrientation convertEdgeDoFOrientation( const EdgeDoFOrientation srcOr, const uint_t v0, const uint_t v1, const uint_t v2 ){
   /// one can calculate the missing local index by adding all indices and subtract from 0+1+2+3
   uint_t v3 = 6 - ( v0 + v1 + v2 );
   switch ( srcOr ){
      case EdgeDoFOrientation::X:
         return getEdgeDoFOrientationFromLocalIDs(v0,v1);
      case EdgeDoFOrientation::Y:
         return getEdgeDoFOrientationFromLocalIDs(v0,v2);
      case EdgeDoFOrientation::XY:
         return getEdgeDoFOrientationFromLocalIDs(v1,v2);
      case EdgeDoFOrientation::Z:
         return getEdgeDoFOrientationFromLocalIDs(v0,v3);
      case EdgeDoFOrientation::XZ:
         return getEdgeDoFOrientationFromLocalIDs(v1,v3);
      case EdgeDoFOrientation::YZ:
         return getEdgeDoFOrientationFromLocalIDs(v2,v3);
      case EdgeDoFOrientation::XYZ:
         /// nothing changes here
         return EdgeDoFOrientation::XYZ;
      default:
         WALBERLA_ABORT("wrong orienation")
   }
}

/// \brief Converts the orientation of an edge DoF in a macro-cell to the respective orientation in a neighboring macro-face.
///
/// \param orientationInCell the edgedof orientation from the cell-local point of view
/// \param cellLocalID0 the first index of the face from cell-local point of view (face-local == 0)
/// \param cellLocalID1 the second index of the face from cell-local point of view (face-local == 1)
/// \param cellLocalID2 the third index of the face from cell-local point of view (face-local == 2)
/// \return the edge DoF orientation from the face-local point of view
inline EdgeDoFOrientation convertEdgeDoFOrientationCellToFace( const EdgeDoFOrientation & orientationInCell, const uint_t & cellLocalID0, const uint_t & cellLocalID1, const uint_t & cellLocalID2 )
{
  uint_t v3 = 6 - ( cellLocalID0 + cellLocalID1 + cellLocalID2 );
  std::map< uint_t, uint_t > cellLocalToFaceLocalIDs;
  cellLocalToFaceLocalIDs[cellLocalID0] = 0;
  cellLocalToFaceLocalIDs[cellLocalID1] = 1;
  cellLocalToFaceLocalIDs[cellLocalID2] = 2;
  cellLocalToFaceLocalIDs[v3] = 3;
  switch ( orientationInCell )
  {
    case EdgeDoFOrientation::X:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[0], cellLocalToFaceLocalIDs[1] );
    case EdgeDoFOrientation::Y:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[0], cellLocalToFaceLocalIDs[2] );
    case EdgeDoFOrientation::XY:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[1], cellLocalToFaceLocalIDs[2] );
    case EdgeDoFOrientation::Z:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[0], cellLocalToFaceLocalIDs[3] );
    case EdgeDoFOrientation::XZ:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[1], cellLocalToFaceLocalIDs[3] );
    case EdgeDoFOrientation::YZ:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalToFaceLocalIDs[2], cellLocalToFaceLocalIDs[3] );
    case EdgeDoFOrientation::XYZ:
      /// nothing changes here
      return EdgeDoFOrientation::XYZ;
    default:
    WALBERLA_ABORT( "wrong orienation" )
  }
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

/// Index of a vertex DoF on a macro face (only access to owned DoFs, no ghost layers).
inline uint_t index( const uint_t & level, const uint_t & x, const uint_t & y, const EdgeDoFOrientation & orientation )
{
  switch ( orientation )
  {
    case EdgeDoFOrientation::X:
      return indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y );
    case EdgeDoFOrientation::Y:
      return 2 * levelToFaceSizeAnyEdgeDoF( level ) + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y );
    case EdgeDoFOrientation::XY:
       return levelToFaceSizeAnyEdgeDoF( level ) + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y );
    default:
      WALBERLA_ASSERT( false, "Invalid orientation" );
      return std::numeric_limits< uint_t >::max();
  }
}

/// Index of a vertex DoF on a ghost layer of a macro face.
/// \param neighbor 0 or 1 for the respective cell neighbor
inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & y, const EdgeDoFOrientation & orientation, const uint_t & neighbor )
{
   uint_t ownDoFs = levelinfo::num_microedges_per_face( level );
   uint_t ghostOnParallelFace = levelinfo::num_microedges_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 1);
   uint_t parallelFaceOneEdgeTypeSize = ghostOnParallelFace / 3;

   /// adjust own dofs if the second cell is used
   uint_t offset = 0;
   if( neighbor == 1 ){
      offset += ownDoFs + ghostOnParallelFace + ghostOnParallelFace / 3;
   }

   switch ( orientation )
   {
      case EdgeDoFOrientation::X:
         return ownDoFs + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
      case EdgeDoFOrientation::Y:
         return ownDoFs + 2 * parallelFaceOneEdgeTypeSize + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
      case EdgeDoFOrientation::XY:
         return ownDoFs + parallelFaceOneEdgeTypeSize + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
      case EdgeDoFOrientation::Z:
         return ownDoFs + ghostOnParallelFace + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y ) + offset;
      case EdgeDoFOrientation::XZ:
         return ownDoFs + ghostOnParallelFace +  2 * levelToFaceSizeAnyEdgeDoF( level ) + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) , x, y ) + offset;
      case EdgeDoFOrientation::YZ:
         return ownDoFs + ghostOnParallelFace + levelToFaceSizeAnyEdgeDoF( level ) + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) , x, y ) + offset;
      case EdgeDoFOrientation::XYZ:
         return ownDoFs * 2 + ghostOnParallelFace + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
      case EdgeDoFOrientation::INVALID:
         return std::numeric_limits< uint_t >::max();
   }
}

inline uint_t horizontalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
   return index( level, col, row, EdgeDoFOrientation::X);
};

inline uint_t verticalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
   return index( level, col, row, EdgeDoFOrientation::Y);
}

inline uint_t diagonalIndex( const uint_t & level, const uint_t & col, const uint_t & row )
{
   return index( level, col, row, EdgeDoFOrientation::XY);
}

// Stencil access functions

inline uint_t indexFromHorizontalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
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

constexpr std::array<stencilDirection ,4> neighborsFromHorizontalEdgeWithoutCenter =
  {{ sD::EDGE_DI_S, sD::EDGE_VE_SE,
     sD::EDGE_DI_N, sD::EDGE_VE_NW
   }};

inline uint_t indexFromDiagonalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
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

constexpr std::array<stencilDirection ,4> neighborsFromDiagonalEdgeWithoutCenter =
  {{ sD::EDGE_HO_S, sD::EDGE_VE_E,
     sD::EDGE_HO_N, sD::EDGE_VE_W
   }};

inline uint_t indexFromVerticalEdge( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
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

constexpr std::array<stencilDirection ,4> neighborsFromVerticalEdgeWithoutCenter =
  {{ sD::EDGE_HO_SE, sD::EDGE_DI_E,
     sD::EDGE_HO_NW, sD::EDGE_DI_W
   }};

inline uint_t indexFromVertex( const uint_t & level, const uint_t & col, const uint_t & row, const stencilDirection & dir )
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
    // assert( false );
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


// ##################
// ### Macro Cell ###
// ##################

namespace macrocell {

/// Direct access functions

typedef stencilDirection sD;

inline constexpr uint_t xIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t yIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
    + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t zIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return 2 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
         + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xyIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return 3 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
         + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xzIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return 4 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
         + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t yzIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return 5 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
         + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xyzIndex( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z )
{
  return 6 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) )
         + indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ) - 1, x, y, z );
}

inline constexpr uint_t index( const uint_t & level, const uint_t & x, const uint_t & y, const uint_t & z, const EdgeDoFOrientation & orientation )
{
  switch ( orientation )
  {
  case EdgeDoFOrientation::X:
    return xIndex( level, x, y, z );
  case EdgeDoFOrientation::Y:
    return yIndex( level, x, y, z );
  case EdgeDoFOrientation::Z:
    return zIndex( level, x, y, z );
  case EdgeDoFOrientation::XY:
    return xyIndex( level, x, y, z );
  case EdgeDoFOrientation::XZ:
    return xzIndex( level, x, y, z );
  case EdgeDoFOrientation::YZ:
    return yzIndex( level, x, y, z );
  case EdgeDoFOrientation::XYZ:
    return xyzIndex( level, x, y, z );
  case EdgeDoFOrientation::INVALID:
    return std::numeric_limits< uint_t >::max();
  }
}


inline bool isInnerXEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 0 ) == 0 && onCellFaces.count( 1 ) == 0;
}

inline bool isInnerYEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 0 ) == 0 && onCellFaces.count( 2 ) == 0;
}

inline bool isInnerZEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 1 ) == 0 && onCellFaces.count( 2 ) == 0;
}

inline bool isInnerXYEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 0 ) == 0 && onCellFaces.count( 3 ) == 0;
}

inline bool isInnerXZEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 1 ) == 0 && onCellFaces.count( 3 ) == 0;
}

inline bool isInnerYZEdgeDoF( const uint_t & level, const indexing::Index & idx )
{
  const auto onCellFaces = indexing::isOnCellFace( idx, levelinfo::num_microvertices_per_edge( level ) - 1 );
  return onCellFaces.count( 2 ) == 0 && onCellFaces.count( 3 ) == 0;
}

inline bool isInnerXYZEdgeDoF( const uint_t &, const indexing::Index & )
{
  return true;
}

inline bool isInnerEdgeDoF( const uint_t & level, const indexing::Index & idx, const EdgeDoFOrientation & orientation )
{
  switch ( orientation )
  {
    case EdgeDoFOrientation::X:
      return isInnerXEdgeDoF( level, idx );
    case EdgeDoFOrientation::Y:
      return isInnerYEdgeDoF( level, idx );
    case EdgeDoFOrientation::Z:
      return isInnerZEdgeDoF( level, idx );
    case EdgeDoFOrientation::XY:
      return isInnerXYEdgeDoF( level, idx );
    case EdgeDoFOrientation::XZ:
      return isInnerXZEdgeDoF( level, idx );
    case EdgeDoFOrientation::YZ:
      return isInnerYZEdgeDoF( level, idx );
    case EdgeDoFOrientation::XYZ:
      return isInnerXYZEdgeDoF( level, idx );
    default:
      WALBERLA_ASSERT( false, "Invalid orientation." );
      return true;
  }
}


// Iterators

class Iterator : public indexing::CellIterator
{
public:
    Iterator( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
      CellIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter )
    {}
};

class IteratorXYZ : public indexing::CellIterator
{
public:
    IteratorXYZ( const uint_t & level, const uint_t & offsetToCenter = 0 ) :
      CellIterator( levelinfo::num_microedges_per_edge( level ) - 1, offsetToCenter )
    {}
};

class BoundaryIterator : public hhg::indexing::CellBorderIterator
{
public:
  BoundaryIterator( const uint_t & level, const uint_t & vertex0, const uint_t & vertex1,
                    const uint_t & vertex2, const uint_t & offsetToCenter = 0 ) :
     CellBorderIterator( levelinfo::num_microedges_per_edge( level ), vertex0, vertex1, vertex2, offsetToCenter )
  {}
};

class BoundaryIteratorXYZ : public hhg::indexing::CellBorderIterator
{
public:
  BoundaryIteratorXYZ( const uint_t & level, const uint_t & vertex0, const uint_t & vertex1,
                    const uint_t & vertex2, const uint_t & offsetToCenter = 0 ) :
     CellBorderIterator( levelinfo::num_microedges_per_edge( level ) - 1, vertex0, vertex1, vertex2, offsetToCenter )
  {}
};

} // namespace macrocell


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

inline bool isHorizontalEdgeOnBoundary(const uint_t level, const hhg::indexing::Index& idx){
  /// level is only needed in the diagonal case
  WALBERLA_UNUSED( level );
  return ( idx.row() == 0 );
}

inline bool isVerticalEdgeOnBoundary(const uint_t level, const hhg::indexing::Index& idx){
  /// level is only needed in the diagonal case
  WALBERLA_UNUSED( level );
  return ( idx.col() == 0 );
}

inline bool isDiagonalEdgeOnBoundary(const uint_t level, const hhg::indexing::Index& idx){
  return ( (idx.col() + idx.row()) == (hhg::levelinfo::num_microedges_per_edge( level ) - 1) );
}


} // namespace edgedof
} // namespace hhg
