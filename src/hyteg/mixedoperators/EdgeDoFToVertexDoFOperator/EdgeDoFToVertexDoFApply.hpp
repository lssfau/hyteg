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

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Elements3D.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace EdgeDoFToVertexDoF {

/// map[neighborCellID][leafOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > > MacroVertexStencilMap_T;
/// map[neighborCellID][leafOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > > MacroEdgeStencilMap_T;
/// map[neighborCellID][leafOrientation][indexOffset] = weight
typedef std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > > MacroFaceStencilMap_T;

/// map[leafOrientation][indexOffset] = weight
typedef std::map< edgedof::EdgeDoFOrientation, std::map< indexing::Index, real_t > > MacroCellStencilMap_T;

inline void applyVertex( uint_t                                                     level,
                         Vertex&                                                    vertex,
                         const PrimitiveDataID< StencilMemory< real_t >, Vertex >&  operatorId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstId,
                         UpdateType                                                 update )
{
   real_t* opr_data = vertex.getData( operatorId )->getPointer( level );
   real_t* src      = vertex.getData( srcId )->getPointer( level );
   real_t* dst      = vertex.getData( dstId )->getPointer( level );

   real_t tmp = 0;
   for ( uint_t i = 0; i < vertex.getData( operatorId )->getSize( level ); ++i )
   {
      tmp += src[i] * opr_data[i];
   }

   if ( update == Replace )
   {
      dst[0] = tmp;
   }
   else if ( update == Add )
   {
      dst[0] += tmp;
   }
}

inline void applyVertex3D( const uint_t&                                                                level,
                           const Vertex&                                                                vertex,
                           const PrimitiveStorage&                                                      storage,
                           const PrimitiveDataID< LevelWiseMemory< MacroVertexStencilMap_T >, Vertex >& operatorId,
                           const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                   srcId,
                           const PrimitiveDataID< FunctionMemory< real_t >, Vertex >&                   dstId,
                           UpdateType                                                                   update )
{
   auto    opr_data = vertex.getData( operatorId )->getData( level );
   real_t* src      = vertex.getData( srcId )->getPointer( level );
   real_t* dst      = vertex.getData( dstId )->getPointer( level );

   const auto centerIndexOnVertex = indexing::Index( 0, 0, 0 );

   real_t tmp = real_c( 0 );

   for ( uint_t neighborCellID = 0; neighborCellID < vertex.getNumNeighborCells(); neighborCellID++ )
   {
      const Cell& neighborCell      = *( storage.getCell( vertex.neighborCells().at( neighborCellID ) ) );
      auto        cellLocalVertexID = neighborCell.getLocalVertexID( vertex.getID() );

      const auto basisInCell = algorithms::getMissingIntegersAscending< 1, 4 >( { cellLocalVertexID } );

      const auto centerIndexInCell = indexing::basisConversion(
          centerIndexOnVertex, basisInCell, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );

      for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientationsWithoutXYZ )
      {
         for ( const auto& stencilIt : opr_data[neighborCellID][leafOrientationInCell] )
         {
            const auto stencilOffset = stencilIt.first;
            const auto stencilWeight = stencilIt.second;

            const auto leafIndexInCell = centerIndexInCell + stencilOffset;

            const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
            const auto onCellEdgesSet = edgedof::macrocell::isOnCellEdges( level, leafIndexInCell, leafOrientationInCell );

            uint_t leafArrayIndexOnVertex = std::numeric_limits< uint_t >::max();

            WALBERLA_ASSERT_GREATER( onCellFacesSet.size(), 0 );

            if ( onCellFacesSet.size() == 1 )
            {
               // on macro-face
               WALBERLA_ASSERT_EQUAL( onCellEdgesSet.size(), 0 );
               const auto faceID            = neighborCell.neighborFaces().at( *onCellFacesSet.begin() );
               const auto vertexLocalFaceID = vertex.face_index( faceID );
               leafArrayIndexOnVertex       = vertex.getNumNeighborEdges() + vertexLocalFaceID;
            }
            else
            {
               // on macro-edge
               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), 2 );
               WALBERLA_ASSERT_EQUAL( onCellEdgesSet.size(), 1 );
               const auto edgeID = neighborCell.neighborEdges().at( *onCellEdgesSet.begin() );
               if ( vertex.neighborPrimitiveExists( edgeID ) )
               {
                  const auto vertexLocalEdgeID = vertex.edge_index( edgeID );
                  leafArrayIndexOnVertex       = vertexLocalEdgeID;
               }
               else
               {
                  // The leaf edgedof is located on an opposite macro-edge
                  // which is no direct neighbor of the current macro-vertex.
                  // This should only happen on level 0.
                  WALBERLA_ASSERT_EQUAL( level, 0 );
                  const std::vector< uint_t > onCellFaces( onCellFacesSet.begin(), onCellFacesSet.end() );
                  WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 2 );
                  const auto faceID0 = neighborCell.neighborFaces().at( onCellFaces.at( 0 ) );
                  const auto faceID1 = neighborCell.neighborFaces().at( onCellFaces.at( 1 ) );

                  WALBERLA_ASSERT( !( vertex.neighborPrimitiveExists( faceID0 ) && vertex.neighborPrimitiveExists( faceID1 ) ) );

                  if ( vertex.neighborPrimitiveExists( faceID0 ) )
                  {
                     const auto vertexLocalFaceID = vertex.face_index( faceID0 );
                     leafArrayIndexOnVertex       = vertex.getNumNeighborEdges() + vertexLocalFaceID;
                  }
                  else if ( vertex.neighborPrimitiveExists( faceID1 ) )
                  {
                     const auto vertexLocalFaceID = vertex.face_index( faceID1 );
                     leafArrayIndexOnVertex       = vertex.getNumNeighborEdges() + vertexLocalFaceID;
                  }
               }
            }

            tmp += src[leafArrayIndexOnVertex] * stencilWeight;
         }
      }
   }

   if ( update == Replace )
   {
      dst[0] = tmp;
   }
   else if ( update == Add )
   {
      dst[0] += tmp;
   }
}

inline void applyEdge( const uint_t&                                            Level,
                       Edge&                                                    edge,
                       const PrimitiveDataID< StencilMemory< real_t >, Edge >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId,
                       UpdateType                                               update )
{
   using namespace hyteg::edgedof;
   size_t rowsize = levelinfo::num_microvertices_per_edge( Level );

   real_t* opr_data = edge.getData( operatorId )->getPointer( Level );
   real_t* src      = edge.getData( srcId )->getPointer( Level );
   real_t* dst      = edge.getData( dstId )->getPointer( Level );

   real_t tmp;

   for ( idx_t i = 1; i < rowsize - 1; ++i )
   {
      tmp = 0.0;
      for ( auto k : macroedge::neighborsOnEdgeFromVertex )
      {
         tmp += opr_data[stencilIndexFromVertex( k )] * src[macroedge::indexFromVertex( Level, i, k )];
      }
      for ( auto k : macroedge::neighborsOnSouthFaceFromVertex )
      {
         tmp += opr_data[stencilIndexFromVertex( k )] * src[macroedge::indexFromVertex( Level, i, k )];
      }
      if ( edge.getNumNeighborFaces() == 2 )
      {
         for ( auto k : macroedge::neighborsOnNorthFaceFromVertex )
         {
            tmp += opr_data[stencilIndexFromVertex( k )] * src[macroedge::indexFromVertex( Level, i, k )];
         }
      }
      if ( update == Replace )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] = tmp;
      }
      else if ( update == Add )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] += tmp;
      }
   }
}

inline void applyEdge3D( const uint_t&                                                            level,
                         const Edge&                                                              edge,
                         const PrimitiveStorage&                                                  storage,
                         const PrimitiveDataID< LevelWiseMemory< MacroEdgeStencilMap_T >, Edge >& operatorId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                 srcId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                 dstId,
                         UpdateType                                                               update )
{
   auto    opr_data = edge.getData( operatorId )->getData( level );
   real_t* src      = edge.getData( srcId )->getPointer( level );
   real_t* dst      = edge.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexOnEdge : hyteg::vertexdof::macroedge::Iterator( level, 1 ) )
   {
      real_t tmp = real_c( 0 );

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             { neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
               neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 ) } );

         const auto centerIndexInCell = indexing::basisConversion(
             centerIndexOnEdge, basisInCell, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ) );

         for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientations )
         {
            for ( const auto& stencilIt : opr_data[neighborCellID][leafOrientationInCell] )
            {
               const auto stencilOffset = stencilIt.first;
               const auto stencilWeight = stencilIt.second;

               const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace(
                   leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               const auto leafIndexInCell = centerIndexInCell + stencilOffset;

               const auto leafIndexOnEdge = leafOrientationOnEdge == edgedof::EdgeDoFOrientation::XYZ ?
                                                edgedof::macrocell::getIndexInNeighboringMacroEdgeXYZ(
                                                    leafIndexInCell, neighborCell, cellLocalEdgeID, storage, level ) :
                                                edgedof::macrocell::getIndexInNeighboringMacroEdge(
                                                    leafIndexInCell, neighborCell, cellLocalEdgeID, storage, level );

               const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
               const auto onCellFacesSetOnEdge =
                   edgedof::macrocell::isOnCellFaces( level, leafIndexOnEdge, leafOrientationOnEdge );

               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

               uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

               const auto& cellLocalIDsOfNeighborFaces =
                   indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
               std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
               std::set_intersection( cellLocalIDsOfNeighborFaces.begin(),
                                      cellLocalIDsOfNeighborFaces.end(),
                                      onCellFacesSet.begin(),
                                      onCellFacesSet.end(),
                                      std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

               if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.empty() )
               {
                  // leaf in macro-cell
                  leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborCell(
                      level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces(), leafOrientationOnEdge );
               }
               else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
               {
                  // leaf on macro-face
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) );
                  const auto cellLocalFaceID = *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin();
                  const auto facePrimitiveID = neighborCell.neighborFaces().at( cellLocalFaceID );
                  WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), facePrimitiveID ) !=
                                   edge.neighborFaces().end() );

                  // The leaf orientation on the edge must be X, Y or XY since it is located on a neighboring face.
                  // Therefore we need to know the three spanning vertex IDs and convert the leaf orientation again.
                  const auto& spanningCellLocalVertices = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( cellLocalFaceID );
                  std::array< uint_t, 4 > faceBasisInCell{};
                  if ( spanningCellLocalVertices.count( basisInCell.at( 2 ) ) == 1 )
                  {
                     faceBasisInCell = basisInCell;
                  }
                  else
                  {
                     WALBERLA_ASSERT( spanningCellLocalVertices.count( basisInCell.at( 3 ) ) == 1 );
                     faceBasisInCell    = basisInCell;
                     faceBasisInCell[2] = basisInCell.at( 3 );
                     faceBasisInCell[3] = basisInCell.at( 2 );
                  }

                  const auto leafIndexOnEdgeGhostLayer = indexing::basisConversion(
                      leafIndexInCell, { 0, 1, 2, 3 }, faceBasisInCell, levelinfo::num_microedges_per_edge( level ) );
                  const auto leafOrientationOnEdgeGhostLayer = edgedof::convertEdgeDoFOrientationCellToFace(
                      leafOrientationInCell, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ), faceBasisInCell.at( 2 ) );

                  const auto localFaceIDOnEdge = edge.face_index( facePrimitiveID );
                  leafArrayIndexOnEdge         = edgedof::macroedge::indexOnNeighborFace(
                      level, leafIndexOnEdgeGhostLayer.x(), localFaceIDOnEdge, leafOrientationOnEdgeGhostLayer );
               }
               else
               {
                  // leaf on macro-edge
                  WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 );
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) );
                  WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, edgedof::EdgeDoFOrientation::X );
                  leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x() );
               }

               tmp += src[leafArrayIndexOnEdge] * stencilWeight;
            }
         }
      }

      if ( update == Replace )
      {
         dst[vertexdof::macroedge::index( level, centerIndexOnEdge.x() )] = tmp;
      }
      else if ( update == Add )
      {
         dst[vertexdof::macroedge::index( level, centerIndexOnEdge.x() )] += tmp;
      }
   }
}

inline void applyFace( const uint_t&                                            Level,
                       Face&                                                    face,
                       const PrimitiveDataID< StencilMemory< real_t >, Face >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                       UpdateType                                               update )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   real_t* opr_data = face.getData( operatorId )->getPointer( Level );
   real_t* src      = face.getData( srcId )->getPointer( Level );
   real_t* dst      = face.getData( dstId )->getPointer( Level );

   real_t tmp;

   using namespace edgedof::macroface;

   for ( idx_t i = 1; i < rowsize - 2; ++i )
   {
      for ( idx_t j = 1; j < inner_rowsize - 2; ++j )
      {
         tmp = 0.0;

         for ( auto k : neighborsFromVertex )
         {
            tmp += opr_data[edgedof::stencilIndexFromVertex( k )] * src[indexFromVertex( Level, i, j, k )];
         }

         if ( update == Replace )
         {
            dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
         }
         else if ( update == Add )
         {
            dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
         }
      }
      --inner_rowsize;
   }
}

inline void applyFace3D( const uint_t&                                                            level,
                         Face&                                                                    face,
                         const PrimitiveStorage&                                                  storage,
                         const PrimitiveDataID< LevelWiseMemory< MacroFaceStencilMap_T >, Face >& operatorId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Face >&                 srcId,
                         const PrimitiveDataID< FunctionMemory< real_t >, Face >&                 dstId,
                         UpdateType                                                               update )
{
   auto    opr_data = face.getData( operatorId )->getData( level );
   real_t* src      = face.getData( srcId )->getPointer( level );
   real_t* dst      = face.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexInFace : hyteg::vertexdof::macroface::Iterator( level, 1 ) )
   {
      std::map< uint_t, real_t > contributionPerCell;

      for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
      {
         contributionPerCell[neighborCellID] = real_c( 0 );

         const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
         const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

         const auto centerIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );

         WALBERLA_ASSERT_GREATER( vertexdof::macrocell::isOnCellFace( centerIndexInCell, level ).size(), 0 );

         for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            for ( const auto& stencilIt : opr_data[neighborCellID][leafOrientation] )
            {
               const auto stencilOffset = stencilIt.first;
               const auto stencilWeight = stencilIt.second;

               const auto leafOrientationInFace = edgedof::macrocell::getOrientationInNeighboringMacroFace(
                   leafOrientation, neighborCell, localFaceID, storage );

               const auto leafIndexInCell = centerIndexInCell + stencilOffset;
               const auto leafIndexInFace = leafOrientation == edgedof::EdgeDoFOrientation::XYZ ?
                                                edgedof::macrocell::getIndexInNeighboringMacroFaceXYZ(
                                                    leafIndexInCell, neighborCell, localFaceID, storage, level ) :
                                                edgedof::macrocell::getIndexInNeighboringMacroFace(
                                                    leafIndexInCell, neighborCell, localFaceID, storage, level );

               WALBERLA_ASSERT_LESS_EQUAL( leafIndexInFace.z(), 1 );

               uint_t leafArrayIndexInFace;
               if ( algorithms::contains( edgedof::faceLocalEdgeDoFOrientations, leafOrientationInFace ) &&
                    leafIndexInFace.z() == 0 )
               {
                  leafArrayIndexInFace =
                      edgedof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace );
               }
               else
               {
                  leafArrayIndexInFace = edgedof::macroface::index(
                      level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace, neighborCellID );
               }

               contributionPerCell[neighborCellID] += stencilWeight * src[leafArrayIndexInFace];
            }
         }
      }

      const auto dstIdx = vertexdof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y() );
      if ( update == Replace )
      {
         dst[dstIdx] = contributionPerCell[0];
         if ( contributionPerCell.size() == 2 )
         {
            dst[dstIdx] += contributionPerCell[1];
         }
      }
      else
      {
         dst[dstIdx] += contributionPerCell[0];
         if ( contributionPerCell.size() == 2 )
         {
            dst[dstIdx] += contributionPerCell[1];
         }
         // WALBERLA_LOG_DEVEL_ON_ROOT( "Edge -> Vertex: " << dst[dstIdx] << " = " << contributionPerCell[0] << " + " << contributionPerCell[1] );
      }
   }
}

inline void applyCell( const uint_t&                                                                                Level,
                       Cell&                                                                                        cell,
                       const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >& operatorId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     srcId,
                       const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                                     dstId,
                       UpdateType                                                                                   update )
{
   auto    opr_data = cell.getData( operatorId )->getData( Level );
   real_t* src      = cell.getData( srcId )->getPointer( Level );
   real_t* dst      = cell.getData( dstId )->getPointer( Level );

   real_t tmp;

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      tmp = 0.0;
      for ( const auto& orientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( orientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            const auto srcIdx      = it + neighbor;
            const auto srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), orientation );
            tmp += opr_data[orientation][neighbor] * src[srcArrayIdx];
         }
      }

      const auto dstArrayIdx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );

      if ( update == Replace )
      {
         dst[dstArrayIdx] = tmp;
      }
      else if ( update == Add )
      {
         dst[dstArrayIdx] += tmp;
      }
   }
}

} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
