/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

inline void saveEdgeOperator( const uint_t&                                           Level,
                              const Edge&                                             edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge >& operatorId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& srcId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& dstId,
                              const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   using namespace hyteg::edgedof::macroedge;
   size_t rowsize = levelinfo::num_microedges_per_edge( Level );

   real_t* opr_data = edge.getData( operatorId )->getPointer( Level );
   idx_t*  src      = edge.getData( srcId )->getPointer( Level );
   idx_t*  dst      = edge.getData( dstId )->getPointer( Level );

   idx_t srcInt;
   idx_t dstInt;

   for ( uint_t i = 0; i < rowsize; ++i )
   {
      dstInt = dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )];

      for ( uint_t k = 0; k < neighborsOnEdgeFromHorizontalEdge.size(); ++k )
      {
         srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnEdgeFromHorizontalEdge[k] )];
         mat->addValue( uint_c( dstInt ),
                        uint_c( srcInt ),
                        opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( neighborsOnEdgeFromHorizontalEdge[k] )] );
      }
      for ( uint_t k = 0; k < neighborsOnSouthFaceFromHorizontalEdge.size(); ++k )
      {
         srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnSouthFaceFromHorizontalEdge[k] )];
         mat->addValue( uint_c( dstInt ),
                        uint_c( srcInt ),
                        opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( neighborsOnSouthFaceFromHorizontalEdge[k] )] );
      }
      if ( edge.getNumNeighborFaces() == 2 )
      {
         for ( uint_t k = 0; k < neighborsOnNorthFaceFromHorizontalEdge.size(); ++k )
         {
            srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnNorthFaceFromHorizontalEdge[k] )];
            mat->addValue(
                uint_c( dstInt ),
                uint_c( srcInt ),
                opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( neighborsOnNorthFaceFromHorizontalEdge[k] )] );
         }
      }
   }
}

inline void saveEdgeIdentityOperator( const uint_t&                                           Level,
                                      const Edge&                                             edge,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& dstId,
                                      const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   using namespace hyteg::edgedof::macroedge;
   size_t rowsize = levelinfo::num_microedges_per_edge( Level );

   idx_t* dst = edge.getData( dstId )->getPointer( Level );
   idx_t  dstInt;

   for ( uint_t i = 0; i < rowsize; ++i )
   {
      dstInt = dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )];
      mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
   }
}

inline void saveEdgeOperator3D( const uint_t&                                                                       level,
                                const Edge&                                                                         edge,
                                const PrimitiveStorage&                                                             storage,
                                const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >& operatorId,
                                const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&                             srcId,
                                const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&                             dstId,
                                const std::shared_ptr< SparseMatrixProxy >&                                         mat )
{
   auto opr_data = edge.getData( operatorId )->getData( level );
   auto src      = edge.getData( srcId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      const edgedof::EdgeDoFOrientation edgeCenterOrientation = edgedof::EdgeDoFOrientation::X;

      const auto dstInt = dst[edgedof::macroedge::index( level, centerIndexOnEdge.x() )];

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             { neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
               neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 ) } );

         const auto centerIndexInCell = indexing::basisConversion(
             centerIndexOnEdge, basisInCell, { 0, 1, 2, 3 }, levelinfo::num_microedges_per_edge( level ) );
         const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
             edgeCenterOrientation, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

         for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientations )
         {
            for ( const auto& stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientationInCell] )
            {
               const auto stencilOffset = stencilIt.first;
               const auto stencilWeight = stencilIt.second;

               const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace(
                   leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               const auto leafIndexInCell = centerIndexInCell + stencilOffset;

               const auto leafIndexOnEdge = indexing::basisConversion(
                   leafIndexInCell, { 0, 1, 2, 3 }, basisInCell, levelinfo::num_microedges_per_edge( level ) );

               const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
               const auto onCellFacesSetOnEdge =
                   edgedof::macrocell::isOnCellFaces( level, leafIndexOnEdge, leafOrientationOnEdge );

               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

               uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

               const auto cellLocalIDsOfNeighborFaces =
                   indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
               std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
               std::set_intersection( cellLocalIDsOfNeighborFaces.begin(),
                                      cellLocalIDsOfNeighborFaces.end(),
                                      onCellFacesSet.begin(),
                                      onCellFacesSet.end(),
                                      std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

               if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
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
                  const auto spanningCellLocalVertices = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( cellLocalFaceID );
                  std::array< uint_t, 4 > faceBasisInCell;
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

               const auto srcInt = src[leafArrayIndexOnEdge];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
            }
         }
      }
   }
}

inline void saveFaceOperator( const uint_t&                                           Level,
                              const Face&                                             face,
                              const PrimitiveDataID< StencilMemory< real_t >, Face >& operatorId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Face >& srcId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstId,
                              const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   real_t* opr_data = face.getData( operatorId )->getPointer( Level );
   idx_t*  src      = face.getData( srcId )->getPointer( Level );
   idx_t*  dst      = face.getData( dstId )->getPointer( Level );

   idx_t srcInt;
   idx_t dstInt;

   using namespace edgedof::macroface;

   for ( const auto& it : hyteg::edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.row() != 0 )
      {
         dstInt = dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
         for ( uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k )
         {
            srcInt = src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
            mat->addValue( uint_c( dstInt ),
                           uint_c( srcInt ),
                           opr_data[edgedof::stencilIndexFromHorizontalEdge( neighborsFromHorizontalEdge[k] )] );
         }
      }
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dstInt = dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
         for ( uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k )
         {
            srcInt = src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
            mat->addValue( uint_c( dstInt ),
                           uint_c( srcInt ),
                           opr_data[edgedof::stencilIndexFromDiagonalEdge( neighborsFromDiagonalEdge[k] )] );
         }
      }
      if ( it.col() != 0 )
      {
         dstInt = dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
         for ( uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k )
         {
            srcInt = src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
            mat->addValue( uint_c( dstInt ),
                           uint_c( srcInt ),
                           opr_data[edgedof::stencilIndexFromVerticalEdge( neighborsFromVerticalEdge[k] )] );
         }
      }
   }
}

inline void saveFaceIdentityOperator( const uint_t&                                           Level,
                                      const Face&                                             face,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Face >& dstId,
                                      const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   idx_t* dst = face.getData( dstId )->getPointer( Level );

   idx_t dstInt;

   using namespace edgedof::macroface;

   for ( const auto& it : hyteg::edgedof::macroface::Iterator( Level, 0 ) )
   {
      if ( it.row() != 0 )
      {
         dstInt = dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
         mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
      }
      if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
         dstInt = dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
         mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
      }
      if ( it.col() != 0 )
      {
         dstInt = dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
         mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
      }
   }
}

inline void saveFaceOperator3D( const uint_t&                                                                       level,
                                const Face&                                                                         face,
                                const PrimitiveStorage&                                                             storage,
                                const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >& operatorId,
                                const PrimitiveDataID< FunctionMemory< idx_t >, Face >&                             srcId,
                                const PrimitiveDataID< FunctionMemory< idx_t >, Face >&                             dstId,
                                const std::shared_ptr< SparseMatrixProxy >&                                         mat )
{
   auto opr_data = face.getData( operatorId )->getData( level );
   auto src      = face.getData( srcId )->getPointer( level );
   auto dst      = face.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexInFace : hyteg::edgedof::macroface::Iterator( level, 0 ) )
   {
      for ( const auto& faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X &&
              edgedof::macroface::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y &&
              edgedof::macroface::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY &&
              edgedof::macroface::isDiagonalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;

         const auto dstIdx =
             edgedof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y(), faceCenterOrientation );
         const auto dstInt = dst[dstIdx];

         for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
         {
            const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
            const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

            const auto centerIndexInCell =
                edgedof::macroface::getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );
            const auto cellCenterOrientation =
                edgedof::macroface::getOrientattionInNeighboringMacroCell( faceCenterOrientation, face, neighborCellID, storage );

            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               for ( const auto& stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientation] )
               {
                  const auto stencilOffset = stencilIt.first;
                  const auto stencilWeight = stencilIt.second;

                  const auto leafOrientationInFace = edgedof::macrocell::getOrientattionInNeighboringMacroFace(
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

                  const auto srcInt = src[leafArrayIndexInFace];
                  mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
               }
            }
         }
      }
   }
}

inline void saveCellOperator( const uint_t&                                                                       Level,
                              const Cell&                                                                         cell,
                              const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >& operatorId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&                             srcId,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Cell >&                             dstId,
                              const std::shared_ptr< SparseMatrixProxy >&                                         mat )
{
   auto srcData  = cell.getData( srcId )->getPointer( Level );
   auto dstData  = cell.getData( dstId )->getPointer( Level );
   auto opr_data = cell.getData( operatorId )->getData( Level );

   for ( const auto& it : edgedof::macrocell::Iterator( Level, 0 ) )
   {
      std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

      if ( edgedof::macrocell::isInnerXEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
      if ( edgedof::macrocell::isInnerYEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
      if ( edgedof::macrocell::isInnerZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
      if ( edgedof::macrocell::isInnerXYEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
      if ( edgedof::macrocell::isInnerXZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
      if ( edgedof::macrocell::isInnerYZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

      for ( const auto& centerOrientation : innerOrientations )
      {
         const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
         const auto dstInt      = dstData[dstArrayIdx];

         for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            const auto edgeDoFNeighbors =
                P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
            for ( const auto& neighbor : edgeDoFNeighbors )
            {
               const auto   srcIdx      = it + neighbor;
               const auto   srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
               const auto   srcInt      = srcData[srcArrayIdx];
               const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
            }
         }
      }
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
   {
      const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
      const auto dstInt      = dstData[dstArrayIdx];

      for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors =
             P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            const auto   srcIdx        = it + neighbor;
            const auto   srcArrayIdx   = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
            const auto   srcInt        = srcData[srcArrayIdx];
            const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
            mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
         }
      }
   }
}

inline void saveCellIdentityOperator( const uint_t&                                           Level,
                                      const Cell&                                             cell,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Cell >& dstId,
                                      const std::shared_ptr< SparseMatrixProxy >&             mat )
{
   auto dstData = cell.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macrocell::Iterator( Level, 0 ) )
   {
      std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

      if ( edgedof::macrocell::isInnerXEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
      if ( edgedof::macrocell::isInnerYEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
      if ( edgedof::macrocell::isInnerZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
      if ( edgedof::macrocell::isInnerXYEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
      if ( edgedof::macrocell::isInnerXZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
      if ( edgedof::macrocell::isInnerYZEdgeDoF( Level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

      for ( const auto& centerOrientation : innerOrientations )
      {
         const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
         const auto dstInt      = dstData[dstArrayIdx];
         mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
      }
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
   {
      const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;

      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );
      const auto dstInt      = dstData[dstArrayIdx];
      mat->addValue( uint_c( dstInt ), uint_c( dstInt ), 1.0 );
   }
}

} // namespace hyteg
