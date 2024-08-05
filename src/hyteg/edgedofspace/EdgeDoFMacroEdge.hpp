/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <algorithm>
#include <vector>

#include "core/math/KahanSummation.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg::edgedof::macroedge {

using indexing::Index;
using walberla::real_c;
using walberla::uint_t;

inline Point3D coordinateFromIndex( const uint_t& level, const Edge& edge, const Index& index )
{
   const real_t  stepFrequency = 1.0 / real_c( levelinfo::num_microedges_per_edge( level ) );
   const Point3D step          = ( edge.getCoordinates()[1] - edge.getCoordinates()[0] ) * stepFrequency;
   return edge.getCoordinates()[0] + 0.5 * step + step * real_c( index.x() );
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               Level,
                         Edge&                                                       edge,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& edgeMemoryId,
                         const ValueType&                                            constant )
{
   auto edgeData = edge.getData( edgeMemoryId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C )] = constant;
   }
}

template < typename ValueType >
inline void interpolate( const std::shared_ptr< PrimitiveStorage >&                                                  storage,
                         const uint_t&                                                                               Level,
                         Edge&                                                                                       edge,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                                 edgeMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >&                  srcIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& ) >& expr )
{
   auto edgeData = edge.getData( edgeMemoryId )->getPointer( Level );

   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( edge.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcIds.size() );

   Point3D xBlend;

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const Point3D currentCoordinates =
          micromesh::microEdgeCenterPosition( storage, edge.getID(), Level, it, edgedof::EdgeDoFOrientation::X );

      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][edgedof::macroedge::horizontalIndex( Level, it.x() )];
      }

      edge.getGeometryMap()->evalF( currentCoordinates, xBlend );
      edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C )] =
          expr( xBlend, srcVector );
   }
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Edge&                                                       edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstID )
{
   auto srcData = edge.getData( srcID );
   auto dstData = edge.getData( dstID );
   srcData->swap( *dstData, level );
}

template < typename ValueType >
inline void add( const uint_t&                                                              Level,
                 Edge&                                                                      edge,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" )
   WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" )

   auto dstData = edge.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      auto tmp = static_cast< ValueType >( 0.0 );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );

      for ( uint_t i = 0; i < scalars.size(); i++ )
      {
         const ValueType scalar  = scalars[i];
         const auto      srcData = edge.getData( srcIds[i] )->getPointer( Level );

         tmp += scalar * srcData[idx];
      }

      dstData[idx] += tmp;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               Level,
                 Edge&                                                       edge,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   auto dstData = edge.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      dstData[idx] += scalar;
   }
}

template < typename ValueType >
inline void assign( const uint_t&                                                              Level,
                    Edge&                                                                      edge,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" )
   WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" )

   auto dstData = edge.getData( dstId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      auto tmp = static_cast< ValueType >( 0.0 );

      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );

      for ( uint_t i = 0; i < scalars.size(); i++ )
      {
         const ValueType scalar  = scalars[i];
         const auto      srcData = edge.getData( srcIds[i] )->getPointer( Level );

         tmp += scalar * srcData[idx];
      }

      dstData[idx] = tmp;
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             Edge&                                                                      edge,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   auto dst = edge.getData( dstId )->getPointer( level );

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );
      ValueType    tmp = edge.getData( srcIds[0] )->getPointer( level )[idx];

      for ( uint_t i = 1; i < srcIds.size(); ++i )
      {
         tmp *= edge.getData( srcIds[i] )->getPointer( level )[idx];
      }

      dst[idx] = tmp;
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               Level,
                      Edge&                                                       edge,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& lhsId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId )
{
   auto lhsData = edge.getData( lhsId )->getPointer( Level );
   auto rhsData = edge.getData( rhsId )->getPointer( Level );

   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      scalarProduct += lhsData[idx] * rhsData[idx];
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                               Level,
                      Edge&                                                       edge,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dataId,
                      const bool&                                                 absolute )
{
   auto data = edge.getData( dataId )->getPointer( Level );

   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      if ( absolute )
      {
         scalarProduct += std::abs( data[idx] );
      }
      else
      {
         scalarProduct += data[idx];
      }
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline void enumerate( const uint_t&                                               Level,
                       Edge&                                                       edge,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                       ValueType&                                                  num )
{
   ValueType* dst = edge.getData( dstId )->getPointer( Level );

   for ( idx_t i = 0; i < levelinfo::num_microedges_per_edge( Level ); ++i )
   {
      dst[hyteg::edgedof::macroedge::horizontalIndex( Level, i )] = num;
      ++num;
   }
}

inline void apply( const uint_t&                                            Level,
                   Edge&                                                    edge,
                   const PrimitiveDataID< StencilMemory< real_t >, Edge >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcId,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId,
                   UpdateType                                               update )
{
   using namespace hyteg::edgedof::macroedge;
   size_t rowsize = levelinfo::num_microedges_per_edge( Level );

   real_t* opr_data = edge.getData( operatorId )->getPointer( Level );
   real_t* src      = edge.getData( srcId )->getPointer( Level );
   real_t* dst      = edge.getData( dstId )->getPointer( Level );

   real_t tmp;

   for ( idx_t i = 0; i < rowsize; ++i )
   {
      tmp = 0.0;
      for ( auto k : neighborsOnEdgeFromHorizontalEdge )
      {
         tmp += opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( k )] * src[indexFromHorizontalEdge( Level, i, k )];
      }
      for ( auto k : neighborsOnSouthFaceFromHorizontalEdge )
      {
         tmp += opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( k )] * src[indexFromHorizontalEdge( Level, i, k )];
      }
      if ( edge.getNumNeighborFaces() == 2 )
      {
         for ( auto k : neighborsOnNorthFaceFromHorizontalEdge )
         {
            tmp += opr_data[hyteg::edgedof::stencilIndexFromHorizontalEdge( k )] * src[indexFromHorizontalEdge( Level, i, k )];
         }
      }

      if ( update == Replace )
      {
         dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] = tmp;
      }
      else if ( update == Add )
      {
         dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] += tmp;
      }
   }
}

inline void apply3D( const uint_t&                                                   level,
                     const Edge&                                                     edge,
                     const PrimitiveStorage&                                         storage,
                     const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Edge >& operatorId,
                     const PrimitiveDataID< FunctionMemory< real_t >, Edge >&        srcId,
                     const PrimitiveDataID< FunctionMemory< real_t >, Edge >&        dstId,
                     UpdateType                                                      update )
{
   auto    opr_data = edge.getData( operatorId )->getData( level );
   real_t* src      = edge.getData( srcId )->getPointer( level );
   real_t* dst      = edge.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      const EdgeDoFOrientation edgeCenterOrientation = EdgeDoFOrientation::X;

      real_t tmp = real_c( 0 );

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

               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() )

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
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) )

                  const auto cellLocalFaceID = *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin();
                  const auto facePrimitiveID = neighborCell.neighborFaces().at( cellLocalFaceID );
                  WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), facePrimitiveID ) !=
                                   edge.neighborFaces().end() )

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
                  WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 )
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) )
                  WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, EdgeDoFOrientation::X )
                  leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x() );
               }

               tmp += src[leafArrayIndexOnEdge] * stencilWeight;
            }
         }
      }

      if ( update == Replace )
      {
         dst[edgedof::macroedge::index( level, centerIndexOnEdge.x() )] = tmp;
      }
      else if ( update == Add )
      {
         dst[edgedof::macroedge::index( level, centerIndexOnEdge.x() )] += tmp;
      }
   }
}

template < typename ValueType >
inline void printFunctionMemory( const uint_t&                                               Level,
                                 const Edge&                                                 edge,
                                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   ValueType* edgeMemory = edge.getData( dstId )->getPointer( Level );
   using namespace std;
   cout << setfill( '=' ) << setw( 100 ) << "" << endl;
   cout << edge << std::left << setprecision( 1 ) << fixed << setfill( ' ' ) << endl;
   uint_t rowsize = levelinfo::num_microvertices_per_edge( Level );
   cout << "Horizontal Edge" << endl;
   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( idx_t i = 1; i < rowsize - 1; ++i )
      {
         cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_NW )]
              << "|";
      }
      cout << endl;
   }
   for ( idx_t i = 1; i < rowsize; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_W )] << "|";
   }
   cout << endl << "     |";
   for ( idx_t i = 1; i < rowsize - 1; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_SE )]
           << "|";
   }
   cout << endl << "Diagonal Edge" << endl;
   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( idx_t i = 1; i < rowsize; ++i )
      {
         cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_DI_NW )]
              << "|";
      }
      cout << endl;
   }
   for ( idx_t i = 0; i < rowsize - 1; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_DI_SE )]
           << "|";
   }
   cout << endl << "Vertical Edge" << endl;
   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( idx_t i = 0; i < rowsize - 1; ++i )
      {
         cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_VE_N )]
              << "|";
      }
      cout << endl;
   }
   for ( idx_t i = 1; i < rowsize; ++i )
   {
      cout << setw( 5 ) << edgeMemory[hyteg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_VE_S )] << "|";
   }
   cout << endl << setfill( '=' ) << setw( 100 ) << "" << endl << setfill( ' ' );
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );
      localMax         = std::max( localMax, src[idx] );
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );
      localMin         = std::min( localMin, src[idx] );
   }

   return localMin;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId )
{
   auto src      = edge.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for ( const auto& it : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C );
      localMax         = std::max( localMax, std::abs( src[idx] ) );
   }

   return localMax;
}

inline void
    invertElementwise( const uint_t& level, Edge& edge, const PrimitiveDataID< FunctionMemory< real_t >, Edge >& edgeDataID )
{
   real_t* data = edge.getData( edgeDataID )->getPointer( level );
   for ( const auto& iter : edgedof::macroedge::Iterator( level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( level, iter.x(), stencilDirection::EDGE_HO_C );
      data[idx]        = real_c( 1.0 ) / data[idx];
   }
}

template < typename ValueType >
inline void createVectorFromFunction( const uint_t&                                               Level,
                                      Edge&                                                       edge,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = edge.getData( srcId )->getPointer( Level );
   auto numerator = edge.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      vec->setValue( numerator[idx], src[idx] );
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               Level,
                                      Edge&                                                       edge,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                                      const PrimitiveDataID< FunctionMemory< idx_t >, Edge >&     numeratorId,
                                      const std::shared_ptr< VectorProxy >&                       vec )
{
   auto src       = edge.getData( srcId )->getPointer( Level );
   auto numerator = edge.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      src[idx]         = vec->getValue( numerator[idx] );
   }
}

inline void applyDirichletBC( const uint_t&                                           Level,
                              Edge&                                                   edge,
                              std::vector< idx_t >&                                   mat,
                              const PrimitiveDataID< FunctionMemory< idx_t >, Edge >& numeratorId )
{
   auto numerator = edge.getData( numeratorId )->getPointer( Level );

   for ( const auto& it : edgedof::macroedge::Iterator( Level ) )
   {
      const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.x(), stencilDirection::EDGE_HO_C );
      mat.push_back( static_cast< idx_t >( numerator[idx] ) );
   }
}

} // namespace hyteg::edgedof::macroedge
