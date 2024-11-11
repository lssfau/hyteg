/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cassert>

#include "hyteg/Levelinfo.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg {
namespace edgedof {

using walberla::uint_c;

constexpr uint_t levelToWidthAnyEdgeDoF( const uint_t& level )
{
   return levelinfo::num_microedges_per_edge( level );
}

constexpr uint_t levelToFaceSizeAnyEdgeDoF( const uint_t& level )
{
   return levelinfo::num_microedges_per_face( level ) / 3;
}

const std::array< EdgeDoFOrientation, 7 > allEdgeDoFOrientations = {
    EdgeDoFOrientation::X,
    EdgeDoFOrientation::Y,
    EdgeDoFOrientation::Z,
    EdgeDoFOrientation::XY,
    EdgeDoFOrientation::XZ,
    EdgeDoFOrientation::YZ,
    EdgeDoFOrientation::XYZ,
};

const std::array< EdgeDoFOrientation, 3 > faceLocalEdgeDoFOrientations = { EdgeDoFOrientation::X,
                                                                           EdgeDoFOrientation::Y,
                                                                           EdgeDoFOrientation::XY };

const std::array< EdgeDoFOrientation, 6 > allEdgeDoFOrientationsWithoutXYZ = { EdgeDoFOrientation::X,
                                                                               EdgeDoFOrientation::Y,
                                                                               EdgeDoFOrientation::Z,
                                                                               EdgeDoFOrientation::XY,
                                                                               EdgeDoFOrientation::XZ,
                                                                               EdgeDoFOrientation::YZ };

const std::map< EdgeDoFOrientation, std::string > edgeDoFOrientationToString = {
    { EdgeDoFOrientation::X, "X" },
    { EdgeDoFOrientation::Y, "Y" },
    { EdgeDoFOrientation::Z, "Z" },
    { EdgeDoFOrientation::XY, "XY" },
    { EdgeDoFOrientation::XZ, "XZ" },
    { EdgeDoFOrientation::YZ, "YZ" },
    { EdgeDoFOrientation::XYZ, "XYZ" },
    { EdgeDoFOrientation::INVALID, "INVALID" },
};

inline std::ostream& operator<<( std::ostream& out, const EdgeDoFOrientation ornt )
{
   const std::string& orntString = edgeDoFOrientationToString.at( ornt );
   return out << orntString;
}

/// \brief Given two logical vertexdof indices, this function returns the appropriate edgedof orientation.
inline EdgeDoFOrientation calcEdgeDoFOrientation( const indexing::Index& vertexIndex0, const indexing::Index& vertexIndex1 )
{
   const indexing::Index offset = vertexIndex1 - vertexIndex0;
   const int             x      = offset.x();
   const int             y      = offset.y();
   const int             z      = offset.z();
   //
   //  WALBERLA_ASSERT_GREATER( x + y + z, 0 );
   //  WALBERLA_ASSERT_LESS_EQUAL( x, 1 );
   //  WALBERLA_ASSERT_LESS_EQUAL( y, 1 );
   //  WALBERLA_ASSERT_LESS_EQUAL( z, 1 );

   if ( x == 1 && y == 0 && z == 0 )
      return EdgeDoFOrientation::X;
   if ( x == -1 && y == 0 && z == 0 )
      return EdgeDoFOrientation::X;

   if ( x == 0 && y == 1 && z == 0 )
      return EdgeDoFOrientation::Y;
   if ( x == 0 && y == -1 && z == 0 )
      return EdgeDoFOrientation::Y;

   if ( x == 0 && y == 0 && z == 1 )
      return EdgeDoFOrientation::Z;
   if ( x == 0 && y == 0 && z == -1 )
      return EdgeDoFOrientation::Z;

   if ( x == -1 && y == 1 && z == 0 )
      return EdgeDoFOrientation::XY;
   if ( x == 1 && y == -1 && z == 0 )
      return EdgeDoFOrientation::XY;

   if ( x == 1 && y == 0 && z == -1 )
      return EdgeDoFOrientation::XZ;
   if ( x == -1 && y == 0 && z == 1 )
      return EdgeDoFOrientation::XZ;

   if ( x == 0 && y == 1 && z == -1 )
      return EdgeDoFOrientation::YZ;
   if ( x == 0 && y == -1 && z == 1 )
      return EdgeDoFOrientation::YZ;

   if ( x == 1 && y == -1 && z == 1 )
      return EdgeDoFOrientation::XYZ;
   if ( x == -1 && y == 1 && z == -1 )
      return EdgeDoFOrientation::XYZ;

   WALBERLA_ASSERT( false, "Invalid index offset." );
   return EdgeDoFOrientation::INVALID;
}

/// \brief Given two logical vertexdof indices, this function returns the logical index of the edgedof inbetween.
/// This function also implicitly calculates the orientation.
inline indexing::Index calcEdgeDoFIndex( const indexing::Index& vertexIndex0, const indexing::Index& vertexIndex1 )
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
      return ( vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1 ) + indexing::Index( 0, -1, 0 );
   case EdgeDoFOrientation::XZ:
      return ( vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1 ) + indexing::Index( 0, 0, -1 );
   case EdgeDoFOrientation::YZ:
      return ( vertexIndex0.y() < vertexIndex1.y() ? vertexIndex0 : vertexIndex1 ) + indexing::Index( 0, 0, -1 );
   case EdgeDoFOrientation::XYZ:
      return ( vertexIndex0.x() < vertexIndex1.x() ? vertexIndex0 : vertexIndex1 ) + indexing::Index( 0, -1, 0 );
   default:
      WALBERLA_ASSERT( false, "Invaild index offset." );
      return indexing::Index(
          std::numeric_limits< int >::max(), std::numeric_limits< int >::max(), std::numeric_limits< int >::max() );
   }
}

/// \brief Given an orientation, this function returns the logical index offsets of the two neighboring vertexdofs on that edge.
/// If the resulting offsets are added to an logical edge index, the resulting logical indices are those of the neighboring vertices.
inline std::array< indexing::Index, 2 > calcNeighboringVertexDoFIndices( const edgedof::EdgeDoFOrientation& orientation )
{
   std::array< indexing::Index, 2 > vertexIndices = {
       indexing::Index( std::numeric_limits< int >::max(), std::numeric_limits< int >::max(), std::numeric_limits< int >::max() ),
       indexing::Index(
           std::numeric_limits< int >::max(), std::numeric_limits< int >::max(), std::numeric_limits< int >::max() ) };
   vertexIndices[0] = indexing::Index( 0, 0, 0 );

   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      vertexIndices[1] = indexing::Index( 1, 0, 0 );
      break;
   case EdgeDoFOrientation::Y:
      vertexIndices[1] = indexing::Index( 0, 1, 0 );
      break;
   case EdgeDoFOrientation::Z:
      vertexIndices[1] = indexing::Index( 0, 0, 1 );
      break;
   case EdgeDoFOrientation::XY:
      vertexIndices[0] = indexing::Index( 1, 0, 0 );
      vertexIndices[1] = indexing::Index( 0, 1, 0 );
      break;
   case EdgeDoFOrientation::XZ:
      vertexIndices[0] = indexing::Index( 1, 0, 0 );
      vertexIndices[1] = indexing::Index( 0, 0, 1 );
      break;
   case EdgeDoFOrientation::YZ:
      vertexIndices[0] = indexing::Index( 0, 1, 0 );
      vertexIndices[1] = indexing::Index( 0, 0, 1 );
      break;
   case EdgeDoFOrientation::XYZ:
      vertexIndices[0] = indexing::Index( 0, 1, 0 );
      vertexIndices[1] = indexing::Index( 1, 0, 1 );
      break;
   default:
      WALBERLA_ASSERT( false, "Invaild index offset." );
   }

   return vertexIndices;
}

/// \brief Given two local vertex ids on the tetrahedron, this function returns the orientation of the edge inbetween those vertices.
/// It also works for macro faces since the first face in the tetrahedron is equivalent to the reference face.
inline EdgeDoFOrientation getEdgeDoFOrientationFromLocalIDs( const uint_t& vertexIndex0, const uint_t& vertexIndex1 )
{
   WALBERLA_ASSERT_UNEQUAL( vertexIndex0, vertexIndex1 );
   WALBERLA_ASSERT_LESS_EQUAL( vertexIndex0, 3 );
   WALBERLA_ASSERT_LESS_EQUAL( vertexIndex1, 3 );
   static indexing::Index vertexIndexII0, vertexIndexII1;
   switch ( vertexIndex0 )
   {
   case 0:
      vertexIndexII0 = { 0, 0, 0 };
      break;
   case 1:
      vertexIndexII0 = { 1, 0, 0 };
      break;
   case 2:
      vertexIndexII0 = { 0, 1, 0 };
      break;
   case 3:
      vertexIndexII0 = { 0, 0, 1 };
      break;
   default:
      WALBERLA_ABORT( "Wrong vertex ID" );
   }
   switch ( vertexIndex1 )
   {
   case 0:
      vertexIndexII1 = { 0, 0, 0 };
      break;
   case 1:
      vertexIndexII1 = { 1, 0, 0 };
      break;
   case 2:
      vertexIndexII1 = { 0, 1, 0 };
      break;
   case 3:
      vertexIndexII1 = { 0, 0, 1 };
      break;
   default:
      WALBERLA_ABORT( "Wrong vertex ID" );
   }

   return calcEdgeDoFOrientation( vertexIndexII0, vertexIndexII1 );
}

/// \brief Converts the orientation of an edge DoF in a macro-edge to the respective orientation in a neighboring macro-face.
///
/// \param orientationInEdge the edgedof orientation from the edge-local point of view
/// \param faceLocalID0 the first index of the edge from face-local point of view (face-local == 0)
/// \param faceLocalID1 the second index of the edge from face-local point of view (face-local == 1)
/// \return the edge DoF orientation from the face-local point of view
inline EdgeDoFOrientation convertEdgeDoFOrientationEdgeToFace( const EdgeDoFOrientation& orientationInEdge,
                                                               const uint_t&             faceLocalID0,
                                                               const uint_t&             faceLocalID1 )
{
   /// one can calculate the missing local index by adding all indices and subtract from 0+1+2
   uint_t faceLocalID2 = 3 - ( faceLocalID0 + faceLocalID1 );
   switch ( orientationInEdge )
   {
   case EdgeDoFOrientation::X:
      return getEdgeDoFOrientationFromLocalIDs( faceLocalID0, faceLocalID1 );
   case EdgeDoFOrientation::Y:
      return getEdgeDoFOrientationFromLocalIDs( faceLocalID0, faceLocalID2 );
   case EdgeDoFOrientation::XY:
      return getEdgeDoFOrientationFromLocalIDs( faceLocalID1, faceLocalID2 );
   default:
      WALBERLA_ABORT( "wrong orienation" )
   }
}

/// \brief converts the edge orientation of the reference face into the correct orientation on the tetrahedron according to the position
/// and the rotation of the face.
inline EdgeDoFOrientation convertEdgeDoFOrientationFaceToCell( const EdgeDoFOrientation srcOr,
                                                               const uint_t             cellLocalID0,
                                                               const uint_t             cellLocalID1,
                                                               const uint_t             cellLocalID2 )
{
   /// one can calculate the missing local index by adding all indices and subtract from 0+1+2+3
   uint_t v3 = 6 - ( cellLocalID0 + cellLocalID1 + cellLocalID2 );
   switch ( srcOr )
   {
   case EdgeDoFOrientation::X:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID0, cellLocalID1 );
   case EdgeDoFOrientation::Y:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID0, cellLocalID2 );
   case EdgeDoFOrientation::XY:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID1, cellLocalID2 );
   case EdgeDoFOrientation::Z:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID0, v3 );
   case EdgeDoFOrientation::XZ:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID1, v3 );
   case EdgeDoFOrientation::YZ:
      return getEdgeDoFOrientationFromLocalIDs( cellLocalID2, v3 );
   case EdgeDoFOrientation::XYZ:
      /// nothing changes here
      return EdgeDoFOrientation::XYZ;
   default:
      WALBERLA_ABORT( "wrong orienation" )
   }
}

/// \brief Converts the orientation of an edge DoF in a macro-cell to the respective orientation in a neighboring macro-face.
///
/// \param orientationInCell the edgedof orientation from the cell-local point of view
/// \param cellLocalID0 the first index of the face from cell-local point of view (face-local == 0)
/// \param cellLocalID1 the second index of the face from cell-local point of view (face-local == 1)
/// \param cellLocalID2 the third index of the face from cell-local point of view (face-local == 2)
/// \return the edge DoF orientation from the face-local point of view
inline EdgeDoFOrientation convertEdgeDoFOrientationCellToFace( const EdgeDoFOrientation& orientationInCell,
                                                               const uint_t&             cellLocalID0,
                                                               const uint_t&             cellLocalID1,
                                                               const uint_t&             cellLocalID2 )
{
   uint_t                            v3 = 6 - ( cellLocalID0 + cellLocalID1 + cellLocalID2 );
   static std::map< uint_t, uint_t > cellLocalToFaceLocalIDs;
   cellLocalToFaceLocalIDs[cellLocalID0] = 0;
   cellLocalToFaceLocalIDs[cellLocalID1] = 1;
   cellLocalToFaceLocalIDs[cellLocalID2] = 2;
   cellLocalToFaceLocalIDs[v3]           = 3;
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
inline constexpr uint_t index( const uint_t& level, const idx_t& x )
{
   return ::hyteg::indexing::macroEdgeIndex( levelToWidthAnyEdgeDoF( level ), x );
}

/// Index of an edge DoF on a ghost layer that is located on a neighbor-face of a macro edge.
/// \param level Refinement level
/// \param x Index on the macro-edge
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
/// \param orientation of the desired micro-edge
inline uint_t
    indexOnNeighborFace( const uint_t& level, const idx_t& x, const uint_t& neighbor, const EdgeDoFOrientation& orientation )
{
   WALBERLA_ASSERT_EQUAL( std::count( faceLocalEdgeDoFOrientations.begin(), faceLocalEdgeDoFOrientations.end(), orientation ),
                          1,
                          "Invalid orientation." );

   WALBERLA_DEBUG_SECTION()
   {
      if ( level == 0 )
      {
         WALBERLA_ASSERT( orientation != EdgeDoFOrientation::X );
      }
   }

   const uint_t numHorizontalDoFsOnEdge       = ::hyteg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );
   const uint_t numHorizontalDoFsOnGhostLayer = ::hyteg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) - 1 );
   const uint_t numOtherTypeDoFsOnGhostLayer  = ::hyteg::indexing::macroEdgeSize( levelToWidthAnyEdgeDoF( level ) );

   switch ( orientation )
   {
   case EdgeDoFOrientation::X: {
      const uint_t offset =
          numHorizontalDoFsOnEdge + neighbor * ( numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer );
      return offset + index( level, x );
   }
   case EdgeDoFOrientation::Y: {
      const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer + numOtherTypeDoFsOnGhostLayer +
                            neighbor * ( numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer );
      return offset + uint_c( x );
   }
   case EdgeDoFOrientation::XY: {
      const uint_t offset = numHorizontalDoFsOnEdge + numHorizontalDoFsOnGhostLayer +
                            neighbor * ( numHorizontalDoFsOnGhostLayer + 2 * numOtherTypeDoFsOnGhostLayer );
      return offset + uint_c( x );
   }
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

/// Index of an edge DoF on a ghost layer that is located on a neighbor-cell of a macro edge.
/// the amount of dofs on the ghost layer are:
/// X: num_microvertices_per_edge - 2
/// XY: num_microvertices_per_edge
/// others: num_microvertices_per_edge - 1
/// \param level Refinement level
/// \param x Index on the macro-edge
/// \param neighbor 0 to access the first neighbor's data, 1 to access second neighbor, ...
/// \param numNeighborFaces Number of neighboring faces
/// \param orientation of the desired micro-edge
inline uint_t indexOnNeighborCell( const uint_t&             level,
                                   const idx_t&              x,
                                   const uint_t&             neighbor,
                                   const uint_t&             numNeighborFaces,
                                   const EdgeDoFOrientation& orientation )
{
   if ( level == 0 )
   {
      WALBERLA_ASSERT_EQUAL( orientation, EdgeDoFOrientation::YZ );
      WALBERLA_ASSERT_EQUAL( x, 0 );
      return 1 + 2 * numNeighborFaces + neighbor;
   }

   const uint_t offsetToFirstCellDoF = levelinfo::num_microedges_per_edge( level ) +
                                       numNeighborFaces * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 );
   const uint_t xGhostOffset        = levelinfo::num_microedges_per_edge( level ) - 2;
   const uint_t xzGhostOffset       = levelinfo::num_microedges_per_edge( level );
   const uint_t otherDoFsoffset     = levelinfo::num_microedges_per_edge( level ) - 1;
   const uint_t offsetNeighborCells = neighbor * ( xGhostOffset + xzGhostOffset + otherDoFsoffset * 5 );
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      return offsetToFirstCellDoF + offsetNeighborCells + index( level, x );
   case EdgeDoFOrientation::Y:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 0 * otherDoFsoffset + index( level, x );
   case EdgeDoFOrientation::Z:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 1 * otherDoFsoffset + index( level, x );
   case EdgeDoFOrientation::XY:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 2 * otherDoFsoffset + index( level, x );
   case EdgeDoFOrientation::XZ:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 3 * otherDoFsoffset + index( level, x );
   case EdgeDoFOrientation::XYZ:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 4 * otherDoFsoffset + index( level, x );
   case EdgeDoFOrientation::YZ:
      return offsetToFirstCellDoF + offsetNeighborCells + xGhostOffset + 5 * otherDoFsoffset + index( level, x );

   default:
      WALBERLA_ASSERT( false, "Invalid orientation" );
      return std::numeric_limits< uint_t >::max();
   }
}

inline uint_t horizontalIndex( const uint_t& level, const idx_t& col )
{
   return index( level, col );
};

inline uint_t horizontalIndex( const uint_t& level, const idx_t& col, const uint_t& neighbor )
{
   return indexOnNeighborFace( level, col, neighbor, EdgeDoFOrientation::X );
};

inline uint_t verticalIndex( const uint_t& level, const idx_t& col, const uint_t& neighbor )
{
   return indexOnNeighborFace( level, col, neighbor, EdgeDoFOrientation::Y );
};

inline uint_t diagonalIndex( const uint_t& level, const idx_t& col, const uint_t& neighbor )
{
   return indexOnNeighborFace( level, col, neighbor, EdgeDoFOrientation::XY );
};

// Stencil access functions

inline uint_t indexFromHorizontalEdge( const uint_t& level, const idx_t& col, const stencilDirection& dir )
{
   // first  neighbor == south
   // second neighbor == north

   switch ( dir )
   {
   case sD::EDGE_HO_C:
      return horizontalIndex( level, col );
   case sD::EDGE_DI_N:
      return diagonalIndex( level, col, 1 );
   case sD::EDGE_DI_S:
      return verticalIndex( level, col, 0 );
   case sD::EDGE_VE_NW:
      return verticalIndex( level, col, 1 );
   case sD::EDGE_VE_SE:
      return diagonalIndex( level, col, 0 );
   default:
      // assert( false );
      return std::numeric_limits< uint_t >::max();
   }
}

inline uint_t indexFromVertex( const uint_t& level, const idx_t& col, const stencilDirection& dir )
{
   // first  neighbor == south
   // second neighbor == north

   switch ( dir )
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
      return verticalIndex( level, col - 1, 0 );
   case sD::EDGE_DI_SE:
      return verticalIndex( level, col, 0 );
   case sD::EDGE_DI_NW:
      return diagonalIndex( level, col - 1, 1 );
   case sD::EDGE_DI_NE:
      return diagonalIndex( level, col, 1 );
   case sD::EDGE_VE_N:
      return verticalIndex( level, col, 1 );
   case sD::EDGE_VE_S:
      return diagonalIndex( level, col - 1, 0 );
   case sD::EDGE_VE_NW:
      return verticalIndex( level, col - 1, 1 );
   case sD::EDGE_VE_SE:
      return diagonalIndex( level, col, 0 );
   default:
      // assert( false );
      return std::numeric_limits< uint_t >::max();
   }
}

constexpr std::array< stencilDirection, 2 > neighborsOnEdgeFromVertex      = { { sD::EDGE_HO_E, sD::EDGE_HO_W } };
constexpr std::array< stencilDirection, 5 > neighborsOnSouthFaceFromVertex = {
    { sD::EDGE_DI_SW, sD::EDGE_VE_S, sD::EDGE_HO_SE, sD::EDGE_DI_SE, sD::EDGE_VE_SE } };
constexpr std::array< stencilDirection, 5 > neighborsOnNorthFaceFromVertex = {
    { sD::EDGE_DI_NE, sD::EDGE_VE_N, sD::EDGE_HO_NW, sD::EDGE_DI_NW, sD::EDGE_VE_NW } };

constexpr std::array< stencilDirection, 1 > neighborsOnEdgeFromHorizontalEdge      = { { sD::EDGE_HO_C } };
constexpr std::array< stencilDirection, 2 > neighborsOnSouthFaceFromHorizontalEdge = { { sD::EDGE_DI_S, sD::EDGE_VE_SE } };
constexpr std::array< stencilDirection, 2 > neighborsOnNorthFaceFromHorizontalEdge = { { sD::EDGE_DI_N, sD::EDGE_VE_NW } };

class Iterator : public indexing::EdgeIterator
{
 public:
   Iterator( const uint_t& level, const uint_t& offsetToCenter = 0, const bool& backwards = false )
   : EdgeIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter, backwards )
   {}
};

} // namespace macroedge

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

typedef stencilDirection sD;

/// Index of an edge DoF on a macro face (only access to owned DoFs, no ghost layers).
inline uint_t index( const uint_t& level, const idx_t& x, const idx_t& y, const EdgeDoFOrientation& orientation )
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

/// Index of an edge DoF on a ghost layer of a macro face.
/// EXAMPLE: number of DoFs in each direction for level 2:
/// - X  on face:  start  0; size 10; last 9
/// - Y  on face:  start 10; size 10; last 19
/// - XY on face:  start 20; size 10; last 29
/// - X  on ghost: start 30; size 6;  last 35
/// - Y  on ghost: start 36; size 6;  last 41
/// - XY on ghost: start 42; size 6;  last 47
/// - Z  on ghost: start 48; size 10; last 57
/// - XZ on ghost: start 58; size 10; last 67
/// - YZ on ghost: start 68; size 10; last 77
/// - XYZon ghost: start 78; size 6;  last 83
/// \param level Refinement level
/// \param x x-index on the macro-face
/// \param y y-Index on the macro-face
/// \param orientation of the desired micro-edge
/// \param neighbor 0 or 1 for the respective cell neighbor
inline constexpr uint_t
    index( const uint_t& level, const idx_t& x, const idx_t& y, const EdgeDoFOrientation& orientation, const uint_t& neighbor )
{
   uint_t ownDoFs = levelinfo::num_microedges_per_face( level );
   uint_t ghostOnParallelFace =
       levelinfo::num_microedges_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 1 );
   uint_t parallelFaceOneEdgeTypeSize = ghostOnParallelFace / 3;

   /// adjust own dofs if the second cell is used
   uint_t offset = 0;
   if ( neighbor == 1 )
   {
      offset += ownDoFs + ghostOnParallelFace + ghostOnParallelFace / 3;
   }

   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      return ownDoFs + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
   case EdgeDoFOrientation::Y:
      return ownDoFs + 2 * parallelFaceOneEdgeTypeSize + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) +
             offset;
   case EdgeDoFOrientation::XY:
      return ownDoFs + parallelFaceOneEdgeTypeSize + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) +
             offset;
   case EdgeDoFOrientation::Z:
      return ownDoFs + ghostOnParallelFace + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y ) + offset;
   case EdgeDoFOrientation::XZ:
      return ownDoFs + ghostOnParallelFace + 2 * levelToFaceSizeAnyEdgeDoF( level ) +
             indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y ) + offset;
   case EdgeDoFOrientation::YZ:
      return ownDoFs + ghostOnParallelFace + levelToFaceSizeAnyEdgeDoF( level ) +
             indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ), x, y ) + offset;
   case EdgeDoFOrientation::XYZ:
      return ownDoFs * 2 + ghostOnParallelFace + indexing::macroFaceIndex( levelToWidthAnyEdgeDoF( level ) - 1, x, y ) + offset;
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

inline uint_t horizontalIndex( const uint_t& level, const idx_t& col, const idx_t& row )
{
   return index( level, col, row, EdgeDoFOrientation::X );
};

inline uint_t verticalIndex( const uint_t& level, const idx_t& col, const idx_t& row )
{
   return index( level, col, row, EdgeDoFOrientation::Y );
}

inline uint_t diagonalIndex( const uint_t& level, const idx_t& col, const idx_t& row )
{
   return index( level, col, row, EdgeDoFOrientation::XY );
}

// Stencil access functions

inline uint_t indexFromHorizontalEdge( const uint_t& level, const idx_t& col, const idx_t& row, const stencilDirection& dir )
{
   switch ( dir )
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

constexpr std::array< stencilDirection, 5 > neighborsFromHorizontalEdge = {
    { sD::EDGE_HO_C, sD::EDGE_DI_S, sD::EDGE_VE_SE, sD::EDGE_DI_N, sD::EDGE_VE_NW } };

constexpr std::array< stencilDirection, 4 > neighborsFromHorizontalEdgeWithoutCenter = {
    { sD::EDGE_DI_S, sD::EDGE_VE_SE, sD::EDGE_DI_N, sD::EDGE_VE_NW } };

inline uint_t indexFromDiagonalEdge( const uint_t& level, const idx_t& col, const idx_t& row, const stencilDirection& dir )
{
   switch ( dir )
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

constexpr std::array< stencilDirection, 5 > neighborsFromDiagonalEdge = {
    { sD::EDGE_DI_C, sD::EDGE_HO_S, sD::EDGE_VE_E, sD::EDGE_HO_N, sD::EDGE_VE_W } };

constexpr std::array< stencilDirection, 4 > neighborsFromDiagonalEdgeWithoutCenter = {
    { sD::EDGE_HO_S, sD::EDGE_VE_E, sD::EDGE_HO_N, sD::EDGE_VE_W } };

inline uint_t indexFromVerticalEdge( const uint_t& level, const idx_t& col, const idx_t& row, const stencilDirection& dir )
{
   switch ( dir )
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

constexpr std::array< stencilDirection, 5 > neighborsFromVerticalEdge = {
    { sD::EDGE_VE_C, sD::EDGE_HO_SE, sD::EDGE_DI_E, sD::EDGE_HO_NW, sD::EDGE_DI_W } };

constexpr std::array< stencilDirection, 4 > neighborsFromVerticalEdgeWithoutCenter = {
    { sD::EDGE_HO_SE, sD::EDGE_DI_E, sD::EDGE_HO_NW, sD::EDGE_DI_W } };

inline uint_t indexFromVertex( const uint_t& level, const idx_t& col, const idx_t& row, const stencilDirection& dir )
{
   // first  neighbor == south
   // second neighbor == north

   switch ( dir )
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

constexpr std::array< stencilDirection, 12 > neighborsFromVertex = { { sD::EDGE_VE_S,
                                                                       sD::EDGE_HO_SE,
                                                                       sD::EDGE_DI_SE,
                                                                       sD::EDGE_VE_SE,
                                                                       sD::EDGE_HO_E,
                                                                       sD::EDGE_DI_NE,
                                                                       sD::EDGE_VE_N,
                                                                       sD::EDGE_HO_NW,
                                                                       sD::EDGE_DI_NW,
                                                                       sD::EDGE_VE_NW,
                                                                       sD::EDGE_HO_W,
                                                                       sD::EDGE_DI_SW } };

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromHorizontalEdge( const stencilDirection dir )
{
   switch ( dir )
   {
   case sD::EDGE_DI_S:
      return 0;
   case sD::EDGE_VE_SE:
      return 1;
   case sD::EDGE_DI_N:
      return 2;
   case sD::EDGE_VE_NW:
      return 3;
   default:
      return std::numeric_limits< size_t >::max();
   }
}

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromDiagonalEdge( const stencilDirection dir )
{
   switch ( dir )
   {
   case sD::EDGE_HO_S:
      return 0;
   case sD::EDGE_VE_E:
      return 1;
   case sD::EDGE_HO_N:
      return 2;
   case sD::EDGE_VE_W:
      return 3;
   default:
      return std::numeric_limits< size_t >::max();
   }
}

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are randomly chosen but need to be kept this way
constexpr inline uint_t stencilIndexFromVerticalEdge( const stencilDirection dir )
{
   switch ( dir )
   {
   case sD::EDGE_HO_S:
      return 0;
   case sD::EDGE_VE_E:
      return 1;
   case sD::EDGE_HO_N:
      return 2;
   case sD::EDGE_VE_W:
      return 3;
   default:
      return std::numeric_limits< size_t >::max();
   }
}

// Iterators

class Iterator : public indexing::FaceIterator
{
 public:
   Iterator( const uint_t& level, const uint_t& offsetToCenter = 0 )
   : FaceIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter )
   {}
};

class BoundaryIterator : public indexing::FaceBoundaryIterator
{
 public:
   BoundaryIterator( const uint_t&                          level,
                     const indexing::FaceBoundaryDirection& direction,
                     const uint_t&                          offsetToCenter     = 0,
                     const uint_t&                          offsetFromVertices = 0 )
   : FaceBoundaryIterator( levelinfo::num_microedges_per_edge( level ), direction, offsetToCenter, offsetFromVertices )
   {}
};

inline bool isHorizontalEdgeOnBoundary( const uint_t level, const hyteg::indexing::Index& idx )
{
   /// level is only needed in the diagonal case
   WALBERLA_UNUSED( level );
   return ( idx.y() == 0 );
}

inline bool isVerticalEdgeOnBoundary( const uint_t level, const hyteg::indexing::Index& idx )
{
   /// level is only needed in the diagonal case
   WALBERLA_UNUSED( level );
   return ( idx.x() == 0 );
}

inline bool isDiagonalEdgeOnBoundary( const uint_t level, const hyteg::indexing::Index& idx )
{
   return ( ( idx.x() + idx.y() ) == ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) );
}

inline bool isInnerEdgeDoF( const uint_t& level, const indexing::Index& idx, const EdgeDoFOrientation& orientation )
{
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      return !isHorizontalEdgeOnBoundary( level, idx );
   case EdgeDoFOrientation::Y:
      return !isVerticalEdgeOnBoundary( level, idx );
   case EdgeDoFOrientation::XY:
      return !isDiagonalEdgeOnBoundary( level, idx );
   default:
      WALBERLA_ASSERT( false, "Invalid orientation." );
      return true;
   }
}

} // namespace macroface

// ##################
// ### Macro Cell ###
// ##################

namespace macrocell {

/// Direct access functions

typedef stencilDirection sD;

inline constexpr uint_t xIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t yIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t zIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return 2 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xyIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return 3 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xzIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return 4 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t yzIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return 5 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ), x, y, z );
}

inline constexpr uint_t xyzIndex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z )
{
   return 6 * levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) +
          indexing::macroCellIndex( levelinfo::num_microedges_per_edge( level ) - 1, x, y, z );
}

inline constexpr uint_t
    index( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z, const EdgeDoFOrientation& orientation )
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
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

inline bool isInnerXEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level > 0 && idx.y() > 0 && idx.z() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level );
}

inline bool isInnerYEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level > 0 && idx.x() > 0 && idx.z() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level );
}

inline bool isInnerZEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level > 0 && idx.x() > 0 && idx.y() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level );
}

inline bool isInnerXYEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level >= 2 && idx.z() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level ) - 1;
}

inline bool isInnerXZEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level >= 2 && idx.y() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level ) - 1;
}

inline bool isInnerYZEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level >= 2 && idx.x() > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level ) - 1;
}

inline bool isInnerXYZEdgeDoF( const uint_t& level, const indexing::Index& idx )
{
   return level > 0 && idx.sum() < levelinfo::num_microedges_per_edge( level ) - 1;
}

inline bool isInnerEdgeDoF( const uint_t& level, const indexing::Index& idx, const EdgeDoFOrientation& orientation )
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

inline std::set< uint_t > isOnCellEdges( const uint_t& level, const indexing::Index& idx, const EdgeDoFOrientation& orientation )
{
   auto               onCellEdgesSet = indexing::isOnCellEdge( idx, levelinfo::num_microedges_per_edge( level ) );
   std::set< uint_t > result;
   uint_t             validLocalEdgeIdx = 7;
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      validLocalEdgeIdx = 0;
      break;
   case EdgeDoFOrientation::Y:
      validLocalEdgeIdx = 1;
      break;
   case EdgeDoFOrientation::Z:
      validLocalEdgeIdx = 3;
      break;
   case EdgeDoFOrientation::XY:
      validLocalEdgeIdx = 2;
      break;
   case EdgeDoFOrientation::XZ:
      validLocalEdgeIdx = 4;
      break;
   case EdgeDoFOrientation::YZ:
      validLocalEdgeIdx = 5;
      break;
   case EdgeDoFOrientation::XYZ:
      break;
   default:
      WALBERLA_ASSERT( false, "Invalid orientation." );
      break;
   }
   if ( onCellEdgesSet.count( validLocalEdgeIdx ) == 1 )
   {
      WALBERLA_ASSERT_LESS_EQUAL( validLocalEdgeIdx, 5 );
      result.insert( validLocalEdgeIdx );
   }
   WALBERLA_ASSERT_LESS_EQUAL( result.size(), 1, "Edgedof cannot lie on more than one edge." );
   return result;
}

inline std::set< uint_t > isOnCellFaces( const uint_t& level, const indexing::Index& idx, const EdgeDoFOrientation& orientation )
{
   auto onCellFacesSet = indexing::isOnCellFace( idx, levelinfo::num_microedges_per_edge( level ) );
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
      onCellFacesSet.erase( 2 );
      onCellFacesSet.erase( 3 );
      break;
   case EdgeDoFOrientation::Y:
      onCellFacesSet.erase( 1 );
      onCellFacesSet.erase( 3 );
      break;
   case EdgeDoFOrientation::Z:
      onCellFacesSet.erase( 0 );
      onCellFacesSet.erase( 3 );
      break;
   case EdgeDoFOrientation::XY:
      onCellFacesSet.erase( 1 );
      onCellFacesSet.erase( 2 );
      break;
   case EdgeDoFOrientation::XZ:
      onCellFacesSet.erase( 0 );
      onCellFacesSet.erase( 2 );
      break;
   case EdgeDoFOrientation::YZ:
      onCellFacesSet.erase( 0 );
      onCellFacesSet.erase( 1 );
      break;
   case EdgeDoFOrientation::XYZ:
      onCellFacesSet.clear();
      break;
   default:
      WALBERLA_ASSERT( false, "Invalid orientation." );
      break;
   }
   WALBERLA_ASSERT_NOT_IDENTICAL( onCellFacesSet.size(), 3, "Edgedof cannot lie in corner of the cell (i.e. on three faces)." );
   return onCellFacesSet;
}

/// \brief Given any edgedof orientation this function returns an arbitrary index of an edgedof that lies in the interior of a macro-cell.
inline indexing::Index getInnerIndexByOrientation( const EdgeDoFOrientation& orientation )
{
   switch ( orientation )
   {
   case EdgeDoFOrientation::X:
   case EdgeDoFOrientation::Y:
   case EdgeDoFOrientation::Z:
      return indexing::Index( 1, 1, 1 );
   case EdgeDoFOrientation::XYZ:
      return indexing::Index( 0, 0, 0 );
   case EdgeDoFOrientation::XY:
      return indexing::Index( 0, 0, 1 );
   case EdgeDoFOrientation::XZ:
      return indexing::Index( 0, 1, 0 );
   case EdgeDoFOrientation::YZ:
      return indexing::Index( 1, 0, 0 );
   default:
      WALBERLA_ASSERT( false, "Invalid orientation." );
      return indexing::Index::max();
   }
}

// Iterators

class Iterator : public indexing::CellIterator
{
 public:
   Iterator( const uint_t& level, const uint_t& offsetToCenter = 0 )
   : CellIterator( levelinfo::num_microedges_per_edge( level ), offsetToCenter )
   {}
};

class IteratorXYZ : public indexing::CellIterator
{
 public:
   IteratorXYZ( const uint_t& level, const uint_t& offsetToCenter = 0 )
   : CellIterator( levelinfo::num_microedges_per_edge( level ) - 1, offsetToCenter )
   {}
};

class BoundaryIterator : public hyteg::indexing::CellBoundaryIterator
{
 public:
   BoundaryIterator( const uint_t& level,
                     const uint_t& vertex0,
                     const uint_t& vertex1,
                     const uint_t& vertex2,
                     const uint_t& offsetToCenter = 0 )
   : CellBoundaryIterator( levelinfo::num_microedges_per_edge( level ), vertex0, vertex1, vertex2, offsetToCenter )
   {}
};

class BoundaryIteratorXYZ : public hyteg::indexing::CellBoundaryIterator
{
 public:
   BoundaryIteratorXYZ( const uint_t& level,
                        const uint_t& vertex0,
                        const uint_t& vertex1,
                        const uint_t& vertex2,
                        const uint_t& offsetToCenter = 0 )
   : CellBoundaryIterator( levelinfo::num_microedges_per_edge( level ) - 1, vertex0, vertex1, vertex2, offsetToCenter )
   {}
};

} // namespace macrocell

/// these numbers specify the postion of each stencil entry in the stencil memory array
/// they are chosen such that the edge dofs on the south face from a macro edge are the first seven entries
/// otherwise the returned index would be out of bounds in the stencil memory array
constexpr inline uint_t stencilIndexFromVertex( const stencilDirection dir )
{
   typedef stencilDirection sD;
   switch ( dir )
   {
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
      return std::numeric_limits< size_t >::max();
   }
}

constexpr inline uint_t stencilIndexFromHorizontalEdge( const stencilDirection dir )
{
   typedef stencilDirection sD;
   switch ( dir )
   {
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
      return std::numeric_limits< size_t >::max();
   }
}

constexpr inline uint_t stencilIndexFromDiagonalEdge( const stencilDirection dir )
{
   typedef stencilDirection sD;
   switch ( dir )
   {
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
      return std::numeric_limits< size_t >::max();
   }
}

constexpr inline uint_t stencilIndexFromVerticalEdge( const stencilDirection dir )
{
   typedef stencilDirection sD;
   switch ( dir )
   {
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
      return std::numeric_limits< size_t >::max();
   }
}

/// Return data indices for the edge dofs of a given micro-cell in a macro-cell
inline void getEdgeDoFDataIndicesFromMicroCell( const indexing::Index&   microCellIndex,
                                                const celldof::CellType& cellType,
                                                const uint_t             level,
                                                std::array< uint_t, 6 >& edgeDoFIndices )
{
   // get indices of micro-vertices forming the micro-cell
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCellIndex, cellType );

   // loop over micro-vertex pairs
   std::array< std::pair< uint_t, uint_t >, 6 > pairs = { std::make_pair( 0, 1 ),
                                                          std::make_pair( 0, 2 ),
                                                          std::make_pair( 0, 3 ),
                                                          std::make_pair( 1, 2 ),
                                                          std::make_pair( 1, 3 ),
                                                          std::make_pair( 2, 3 ) };

   for ( uint_t k = 0; k < 6; ++k )
   {
      // generate Index for vertices in the current pair
      uint_t n1 = pairs[k].first;
      uint_t n2 = pairs[k].second;

      indexing::Index vertexIdx1( (int) verts[n1].x(), (int) verts[n1].y(), (int) verts[n1].z() );
      indexing::Index vertexIdx2( (int) verts[n2].x(), (int) verts[n2].y(), (int) verts[n2].z() );

      // get edge orientation
      EdgeDoFOrientation orientation = calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );
      // edgedof::EdgeDoFOrientation orientation = edgedof::calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );

      // get index for this edgeDoF
      indexing::Index eIdx = edgedof::calcEdgeDoFIndex( vertexIdx1, vertexIdx2 );

      // finally get the data index of the edgeDoF
      edgeDoFIndices[k] = macrocell::index( level, eIdx.x(), eIdx.y(), eIdx.z(), orientation );
   }
}

/// Return data indices for the edge dofs of a given micro-face in a macro-face
inline void getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( const indexing::Index&   microFaceIndex,
                                                              const facedof::FaceType& faceType,
                                                              const uint_t             level,
                                                              std::array< uint_t, 3 >& edgeDoFIndices )
{
   // get indices of micro-vertices forming the micro-cell
   std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace( microFaceIndex, faceType );

   // loop over micro-vertex pairs
   std::array< std::pair< uint_t, uint_t >, 3 > pairs = {
       std::make_pair( 1, 2 ),
       std::make_pair( 0, 2 ),
       std::make_pair( 0, 1 ),
   };

   for ( uint_t k = 0; k < 3; ++k )
   {
      // generate Indexs for vertices in the current pair
      uint_t n1 = pairs[k].first;
      uint_t n2 = pairs[k].second;

      indexing::Index vertexIdx1( (int) verts[n1].x(), (int) verts[n1].y(), 0 );
      indexing::Index vertexIdx2( (int) verts[n2].x(), (int) verts[n2].y(), 0 );

      // get edge orientation
      EdgeDoFOrientation orientation = calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );

      // get index for this edgeDoF
      indexing::Index eIdx = edgedof::calcEdgeDoFIndex( vertexIdx1, vertexIdx2 );

      // finally get the data index of the edgeDoF
      edgeDoFIndices[k] = macroface::index( level, eIdx.x(), eIdx.y(), orientation );
   }
}

/// Return data indices for the edge dofs of a given micro-cell in a macro-cell
inline void getEdgeDoFDataIndicesFromMicroVerticesFEniCSOrdering( const std::array< indexing::Index, 4 >& verts,
                                                                  const uint_t                            level,
                                                                  std::array< uint_t, 6 >&                edgeDoFIndices )
{
   // loop over micro-vertex pairs
   std::array< std::pair< uint_t, uint_t >, 6 > pairs = { std::make_pair( 2, 3 ),
                                                          std::make_pair( 1, 3 ),
                                                          std::make_pair( 1, 2 ),
                                                          std::make_pair( 0, 3 ),
                                                          std::make_pair( 0, 2 ),
                                                          std::make_pair( 0, 1 ) };

   for ( uint_t k = 0; k < 6; ++k )
   {
      // generate Indexs for vertices in the current pair
      uint_t n1 = pairs[k].first;
      uint_t n2 = pairs[k].second;

      indexing::Index vertexIdx1( (int) verts[n1].x(), (int) verts[n1].y(), (int) verts[n1].z() );
      indexing::Index vertexIdx2( (int) verts[n2].x(), (int) verts[n2].y(), (int) verts[n2].z() );

      // get edge orientation
      EdgeDoFOrientation orientation = calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );
      // edgedof::EdgeDoFOrientation orientation = edgedof::calcEdgeDoFOrientation( vertexIdx1, vertexIdx2 );

      // get index for this edgeDoF
      indexing::Index eIdx = edgedof::calcEdgeDoFIndex( vertexIdx1, vertexIdx2 );

      // finally get the data index of the edgeDoF
      edgeDoFIndices[k] = macrocell::index( level, eIdx.x(), eIdx.y(), eIdx.z(), orientation );
   }
}

/// Return data indices for the edge dofs of a given micro-cell in a macro-cell
/// \sa hyteg::n1e1::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering
inline void getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( const indexing::Index&   microCellIndex,
                                                              const celldof::CellType& cellType,
                                                              const uint_t             level,
                                                              std::array< uint_t, 6 >& edgeDoFIndices )
{
   // get indices of micro-vertices forming the micro-cell
   const std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCellIndex, cellType );
   getEdgeDoFDataIndicesFromMicroVerticesFEniCSOrdering( verts, level, edgeDoFIndices );
}

} // namespace edgedof
} // namespace hyteg
