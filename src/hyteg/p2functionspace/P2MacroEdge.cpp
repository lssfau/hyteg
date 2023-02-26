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
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroedge_P2_update_vertexdofs.hpp"

#include "P2MacroFace.hpp"

#include "core/OpenMP.h"

namespace hyteg {
namespace P2 {
namespace macroedge {

void smoothSOR( const uint_t&                                            level,
                const Edge&                                              edge,
                const real_t&                                            relax,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToVertexStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToVertexStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  vertexToEdgeStencilID,
                const PrimitiveDataID< StencilMemory< real_t >, Edge >&  edgeToEdgeStencilID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID )
{
   real_t* vertexToVertexStencil = edge.getData( vertexToVertexStencilID )->getPointer( level );
   real_t* edgeToVertexStencil   = edge.getData( edgeToVertexStencilID )->getPointer( level );
   real_t* dstVertexDoF          = edge.getData( dstVertexDoFID )->getPointer( level );
   real_t* vertexToEdgeStencil   = edge.getData( vertexToEdgeStencilID )->getPointer( level );
   real_t* edgeToEdgeStencil     = edge.getData( edgeToEdgeStencilID )->getPointer( level );
   real_t* dstEdgeDoF            = edge.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhsVertexDoF          = edge.getData( rhsVertexDoFID )->getPointer( level );
   real_t* rhsEdgeDoF            = edge.getData( rhsEdgeDoFID )->getPointer( level );

   real_t tmpVertex = 0, tmpEdgeHO = 0;

   const real_t invVertexCenter = real_c( 1.0 ) / vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
   const real_t invEdgeXCenter  = real_c( 1.0 ) / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )];

   for( const auto& it : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      ////////// VERTEX //////////
      if( it.x() != 0 )
      {
         tmpVertex = rhsVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )];
         /// on edge vertex dof
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
         {
            tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// on edge edge dof
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnEdgeFromVertex )
         {
            tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         /// south face vertex dof
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
         {
            tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// south face edge
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnSouthFaceFromVertex )
         {
            tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         if( edge.getNumNeighborFaces() == 2 )
         {
            /// north face vertex dof
            for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
            {
               tmpVertex -= dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                            vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
            }
            /// north face edge
            for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnNorthFaceFromVertex )
            {
               tmpVertex -= dstEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                            edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
            }
         }
         dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )] =
             (real_c( 1.0 ) - relax) * dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )] +
             relax * invVertexCenter * tmpVertex;
      }
      ////////// HORIZONTAL EDGE //////////
      tmpEdgeHO = rhsEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )];
      /// on edge
      for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on south face
      for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }

      for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge )
      {
         tmpEdgeHO -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on north face
      if( edge.getNumNeighborFaces() == 2 )
      {
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
         {
            tmpEdgeHO -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
         }
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge )
         {
            tmpEdgeHO -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
         }
      }
      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )] =
          (real_c( 1.0 ) - relax) * dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )] +
          relax * invEdgeXCenter * tmpEdgeHO;
   }
}

#ifdef HYTEG_USE_GENERATED_KERNELS
static void smoothSOR3DUpdateVertexDoFsGenerated(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    Edge&                                                                                        edge,
    const real_t&                                                                                relax,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge >&        vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards )
{
   using edgedof::EdgeDoFOrientation;
   using indexing::Index;

   auto v2v_operator = edge.getData( vertexToVertexOperatorId )->getData( level );
   auto e2v_operator = edge.getData( edgeToVertexOperatorId )->getData( level );
   auto v2e_operator = edge.getData( vertexToEdgeOperatorId )->getData( level );
   auto e2e_operator = edge.getData( edgeToEdgeOperatorId )->getData( level );

   real_t* vertexDoFDst = edge.getData( vertexDoFDstId )->getPointer( level );
   real_t* vertexDoFRhs = edge.getData( vertexDoFRhsId )->getPointer( level );
   real_t* edgeDoFDst   = edge.getData( edgeDoFDstId )->getPointer( level );
   real_t* edgeDoFRhs   = edge.getData( edgeDoFRhsId )->getPointer( level );

   real_t center = 0;
   for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
   {
      center += v2v_operator[neighborCellID][indexing::Index( {0, 0, 0} )];
   }

   const real_t vertexDoFRelaxOverCenter = relax / center;
   const real_t oneMinusRelax            = real_c( 1 ) - relax;

   WALBERLA_UNUSED( edgeDoFRhs );

   // pre-calculate for each neighbor cell:
   std::vector< std::array< uint_t, 4 > > cellLocalVertexIDs;
   std::vector< uint_t >                  edgeLocalFace0ID;
   std::vector< uint_t >                  edgeLocalFace1ID;
   for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
   {
      const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
      auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

     // 1. the cell local vertex IDs as seen from the edge locally
     auto cellLocalVertexIDsEntry = algorithms::getMissingIntegersAscending< 2, 4 >(
      {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
       neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );
      cellLocalVertexIDs.push_back( cellLocalVertexIDsEntry );

      // 2. the edge local ID of the lower face
      auto cellLocalIDLowerFace = indexing::getCellLocalFaceIDFromCellLocalVertexIDs(
      cellLocalVertexIDsEntry[0], cellLocalVertexIDsEntry[1], cellLocalVertexIDsEntry[2] );
      auto lowerFacePrimitiveID = neighborCell.neighborFaces().at( cellLocalIDLowerFace );
      edgeLocalFace0ID.push_back( edge.face_index( lowerFacePrimitiveID ) );

      // 3. the edge local ID of the upper face
      auto cellLocalIDUpperFace = indexing::getCellLocalFaceIDFromCellLocalVertexIDs(
      cellLocalVertexIDsEntry[0], cellLocalVertexIDsEntry[1], cellLocalVertexIDsEntry[3] );
      auto upperFacePrimitiveID = neighborCell.neighborFaces().at( cellLocalIDUpperFace );
      edgeLocalFace1ID.push_back( edge.face_index( upperFacePrimitiveID ) );
   }

   for ( const auto& centerIndexOnEdge : hyteg::vertexdof::macroedge::Iterator( level, 1, backwards ) )
   {
      const auto dstIdx     = vertexdof::macroedge::index( level, centerIndexOnEdge.x() );
      real_t     stencilSum = vertexDoFRhs[dstIdx];

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         real_t partialStencilSum = 0.0;

         P2::macroedge::generated::sor_3D_macroedge_P2_update_vertexdofs(
             edgeDoFDst,
             &partialStencilSum,
             vertexDoFDst,
             e2v_operator[neighborCellID],
             static_cast< int64_t >( neighborCellID ),
             static_cast< int64_t >( edgeLocalFace0ID[neighborCellID] ),
             static_cast< int64_t >( edgeLocalFace1ID[neighborCellID] ),
             static_cast< int32_t >( level ),
             static_cast< int64_t >( centerIndexOnEdge.x() ),
             static_cast< int64_t >( cellLocalVertexIDs[neighborCellID][0] ),
             static_cast< int64_t >( cellLocalVertexIDs[neighborCellID][1] ),
             static_cast< int64_t >( cellLocalVertexIDs[neighborCellID][2] ),
             static_cast< int64_t >( edge.getNumNeighborFaces() ),
             v2v_operator[neighborCellID] );

         stencilSum -= partialStencilSum;
      }

      vertexDoFDst[dstIdx] = oneMinusRelax * vertexDoFDst[dstIdx] + vertexDoFRelaxOverCenter * stencilSum;
   }
}
#endif

static void smoothSOR3DUpdateVertexDoFs(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    Edge&                                                                                        edge,
    const real_t&                                                                                relax,
    const PrimitiveDataID< StencilMemory< real_t >, Edge >&                                      vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards )
{
  using edgedof::EdgeDoFOrientation;
  using indexing::Index;
  typedef stencilDirection sD;

  auto v2v_operator = edge.getData( vertexToVertexOperatorId )->getPointer( level );
  auto e2v_operator = edge.getData( edgeToVertexOperatorId )->getData( level );
  auto v2e_operator = edge.getData( vertexToEdgeOperatorId )->getData( level );
  auto e2e_operator = edge.getData( edgeToEdgeOperatorId )->getData( level );

  real_t *vertexDoFDst = edge.getData( vertexDoFDstId )->getPointer( level );
  real_t *vertexDoFRhs = edge.getData( vertexDoFRhsId )->getPointer( level );
  real_t *edgeDoFDst = edge.getData( edgeDoFDstId )->getPointer( level );
  real_t *edgeDoFRhs = edge.getData( edgeDoFRhsId )->getPointer( level );

  const real_t vertexDoFRelaxOverCenter = relax / v2v_operator[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
  const real_t oneMinusRelax = real_c( 1 ) - relax;

  real_t tmp;

  WALBERLA_UNUSED( edgeDoFRhs );

  // updating vertex unknowns
  for ( const auto & centerIndexOnEdge : hyteg::vertexdof::macroedge::Iterator( level, 1, backwards ))
  {
    const auto dstIdx = vertexdof::macroedge::index( level, centerIndexOnEdge.x() );
    tmp = vertexDoFRhs[ dstIdx ];

    // vertex leaves
    const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
    const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

    const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, centerIndexOnEdge.x(), sD::VERTEX_W );
    const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, centerIndexOnEdge.x(), sD::VERTEX_E );

    tmp -= v2v_operator[ stencilIdxW ] * vertexDoFDst[ dofIdxW ] + v2v_operator[ stencilIdxE ] * vertexDoFDst[ dofIdxE ];

    for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
    {
      const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
      const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
      const auto stencilWeightW = v2v_operator[ stencilIdxWNeighborFace ];
      const auto stencilWeightE = v2v_operator[ stencilIdxENeighborFace ];
      const auto dofIdxWNeighborFace = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, centerIndexOnEdge.x(), neighborFace, sD::VERTEX_W );
      const auto dofIdxENeighborFace = vertexdof::macroedge::indexFromVertexOnNeighborFace( level, centerIndexOnEdge.x(), neighborFace, sD::VERTEX_E );
      tmp -= stencilWeightW * vertexDoFDst[dofIdxWNeighborFace] + stencilWeightE * vertexDoFDst[dofIdxENeighborFace];
    }

    for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
    {
      const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
      const auto dofIdx = vertexdof::macroedge::indexFromVertexOnNeighborCell( level, centerIndexOnEdge.x(), neighborCell, edge.getNumNeighborFaces() );
      tmp -= v2v_operator[ stencilIdx ] * vertexDoFDst[ dofIdx ];
    }

    // edge leaves
    for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
    {
      const Cell & neighborCell = *( storage.getCell( edge.neighborCells().at( neighborCellID )));
      auto cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID());

      const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >( { neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
                                                                                  neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 ) } );

      const auto centerIndexInCell = indexing::basisConversion( centerIndexOnEdge, basisInCell, { 0, 1, 2, 3 }, levelinfo::num_microvertices_per_edge( level ));

      for ( const auto & leafOrientationInCell : edgedof::allEdgeDoFOrientations )
      {
        for ( const auto & stencilIt : e2v_operator[neighborCellID][leafOrientationInCell] )
        {
          const auto stencilOffset = stencilIt.first;
          const auto stencilWeight = stencilIt.second;

          const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace(
              leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
          const auto leafIndexInCell = centerIndexInCell + stencilOffset;

          const auto leafIndexOnEdge = leafOrientationOnEdge == edgedof::EdgeDoFOrientation::XYZ
                                       ? edgedof::macrocell::getIndexInNeighboringMacroEdgeXYZ( leafIndexInCell, neighborCell, cellLocalEdgeID, storage, level )
                                       : edgedof::macrocell::getIndexInNeighboringMacroEdge( leafIndexInCell, neighborCell, cellLocalEdgeID, storage, level );

          const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
          const auto onCellFacesSetOnEdge = edgedof::macrocell::isOnCellFaces( level, leafIndexOnEdge, leafOrientationOnEdge );

          WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size());

          uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

          const auto cellLocalIDsOfNeighborFaces = indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
          std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
          std::set_intersection( cellLocalIDsOfNeighborFaces.begin(), cellLocalIDsOfNeighborFaces.end(),
                                 onCellFacesSet.begin(), onCellFacesSet.end(), std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ));

          if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
          {
            // leaf in macro-cell
            leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborCell( level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces(), leafOrientationOnEdge );
          } else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
          {
            // leaf on macro-face
            WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ));
            const auto cellLocalFaceID = *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin();
            const auto facePrimitiveID = neighborCell.neighborFaces().at( cellLocalFaceID );
            WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), facePrimitiveID ) != edge.neighborFaces().end());

            // The leaf orientation on the edge must be X, Y or XY since it is located on a neighboring face.
            // Therefore we need to know the three spanning vertex IDs and convert the leaf orientation again.
            const auto spanningCellLocalVertices = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( cellLocalFaceID );
            std::array< uint_t, 4 > faceBasisInCell;
            if ( spanningCellLocalVertices.count( basisInCell.at( 2 )) == 1 )
            {
              faceBasisInCell = basisInCell;
            } else
            {
              WALBERLA_ASSERT( spanningCellLocalVertices.count( basisInCell.at( 3 )) == 1 );
              faceBasisInCell = basisInCell;
              faceBasisInCell[2] = basisInCell.at( 3 );
              faceBasisInCell[3] = basisInCell.at( 2 );
            }

            const auto leafIndexOnEdgeGhostLayer = indexing::basisConversion( leafIndexInCell, { 0, 1, 2, 3 }, faceBasisInCell, levelinfo::num_microedges_per_edge( level ));
            const auto leafOrientationOnEdgeGhostLayer = edgedof::convertEdgeDoFOrientationCellToFace( leafOrientationInCell, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ),
                                                                                                       faceBasisInCell.at( 2 ));

            const auto localFaceIDOnEdge = edge.face_index( facePrimitiveID );
            leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborFace( level, leafIndexOnEdgeGhostLayer.x(), localFaceIDOnEdge, leafOrientationOnEdgeGhostLayer );
          } else
          {
            // leaf on macro-edge
            WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 );
            WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ));
            WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, edgedof::EdgeDoFOrientation::X );
            leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x());
          }

          tmp -= edgeDoFDst[leafArrayIndexOnEdge] * stencilWeight;
        }
      }
    }

    vertexDoFDst[ dstIdx ] = oneMinusRelax * vertexDoFDst[ dstIdx ] + vertexDoFRelaxOverCenter * tmp;
  }
}

static void smoothSOR3DUpdateEdgeDoFs(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    Edge&                                                                                        edge,
    const real_t&                                                                                relax,
    const PrimitiveDataID< StencilMemory< real_t >, Edge >&                                      vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards )
{
   using edgedof::EdgeDoFOrientation;
   using indexing::Index;

   auto v2v_operator = edge.getData( vertexToVertexOperatorId )->getPointer( level );
   auto e2v_operator = edge.getData( edgeToVertexOperatorId )->getData( level );
   auto v2e_operator = edge.getData( vertexToEdgeOperatorId )->getData( level );
   auto e2e_operator = edge.getData( edgeToEdgeOperatorId )->getData( level );

   real_t* vertexDoFDst = edge.getData( vertexDoFDstId )->getPointer( level );
   real_t* vertexDoFRhs = edge.getData( vertexDoFRhsId )->getPointer( level );
   real_t* edgeDoFDst   = edge.getData( edgeDoFDstId )->getPointer( level );
   real_t* edgeDoFRhs   = edge.getData( edgeDoFRhsId )->getPointer( level );

   const real_t vertexDoFRelaxOverCenter = relax / v2v_operator[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
   const real_t oneMinusRelax            = real_c( 1 ) - relax;

   real_t tmp;

   WALBERLA_UNUSED( vertexDoFRhs );
   WALBERLA_UNUSED( vertexDoFRelaxOverCenter );

   // updating edge unknowns
   for ( const auto& centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0, backwards ) )
   {
      const EdgeDoFOrientation edgeCenterOrientation = EdgeDoFOrientation::X;

      const auto dstIdx = edgedof::macroedge::index( level, centerIndexOnEdge.x() );
      tmp               = edgeDoFRhs[dstIdx];

      real_t e2eDiagonalEntry = 0;

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
              neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

         const auto centerIndexInCell = indexing::basisConversion(
             centerIndexOnEdge, basisInCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );
         const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
             edgeCenterOrientation, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

         // vertex leaves
         for ( const auto& stencilIt : v2e_operator[neighborCellID][cellCenterOrientation] )
         {
            const auto stencilOffset = stencilIt.first;
            const auto stencilWeight = stencilIt.second;

            const auto leafIndexInCell = centerIndexInCell + stencilOffset;
            const auto leafIndexOnEdge = indexing::basisConversion(
                leafIndexInCell, { 0, 1, 2, 3 }, basisInCell, levelinfo::num_microvertices_per_edge( level ) );

            const auto onCellFacesSet       = vertexdof::macrocell::isOnCellFace( leafIndexInCell, level );
            const auto onCellFacesSetOnEdge = vertexdof::macrocell::isOnCellFace( leafIndexOnEdge, level );

            WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() );

            const auto cellLocalIDsOfNeighborFaces = indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
            std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
            std::set_intersection( cellLocalIDsOfNeighborFaces.begin(),
                                   cellLocalIDsOfNeighborFaces.end(),
                                   onCellFacesSet.begin(),
                                   onCellFacesSet.end(),
                                   std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

            uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

            if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
            {
               // leaf in macro-cell
               leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborCell(
                   level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces() );
            }
            else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
            {
               // leaf on macro-face

               const auto faceID = neighborCell.neighborFaces().at( *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin() );
               WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), faceID ) !=
                                edge.neighborFaces().end() );

               const auto localFaceIDOnEdge = edge.face_index( faceID );
               leafArrayIndexOnEdge = vertexdof::macroedge::indexOnNeighborFace( level, leafIndexOnEdge.x(), localFaceIDOnEdge );
            }
            else
            {
               // leaf on macro-edge
               WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 );
               leafArrayIndexOnEdge = vertexdof::macroedge::index( level, leafIndexOnEdge.x() );
            }

            tmp -= vertexDoFDst[leafArrayIndexOnEdge] * stencilWeight;
         }

         // edge leaves
         for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientations )
         {
            for ( const auto& stencilIt : e2e_operator[neighborCellID][cellCenterOrientation][leafOrientationInCell] )
            {
               const auto stencilOffset = stencilIt.first;
               const auto stencilWeight = stencilIt.second;

               if ( leafOrientationInCell == cellCenterOrientation && stencilOffset == Index( 0, 0, 0 ) )
               {
                  e2eDiagonalEntry += stencilWeight;
                  continue;
               }

               const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace(
                   leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               const auto leafIndexInCell = centerIndexInCell + stencilOffset;

               const auto leafIndexOnEdge = indexing::basisConversion(
                   leafIndexInCell, {0, 1, 2, 3}, basisInCell, levelinfo::num_microedges_per_edge( level ) );

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
                      leafIndexInCell, {0, 1, 2, 3}, faceBasisInCell, levelinfo::num_microedges_per_edge( level ) );
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
                  WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, EdgeDoFOrientation::X );
                  leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x() );
               }

               tmp -= edgeDoFDst[leafArrayIndexOnEdge] * stencilWeight;
            }
         }
      }

      edgeDoFDst[dstIdx] = oneMinusRelax * edgeDoFDst[dstIdx] + ( relax / e2eDiagonalEntry ) * tmp;
   }
}

void smoothSOR3D(
    const uint_t&                                                                                level,
    const PrimitiveStorage&                                                                      storage,
    Edge&                                                                                        edge,
    const real_t&                                                                                relax,
    const PrimitiveDataID< StencilMemory< real_t >, Edge >&                                      vertexToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge >&        vertexToVertexOperatorMapId,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >& edgeToVertexOperatorId,
    const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge >& vertexToEdgeOperatorId,
    const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&          edgeToEdgeOperatorId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     vertexDoFRhsId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFDstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Edge >&                                     edgeDoFRhsId,
    const bool&                                                                                  backwards )
{
   WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "VertexDoFs" ); }
   if ( globalDefines::useGeneratedKernels )
   {
#ifdef HYTEG_USE_GENERATED_KERNELS
      smoothSOR3DUpdateVertexDoFsGenerated( level,
                                            storage,
                                            edge,
                                            relax,
                                            vertexToVertexOperatorMapId,
                                            edgeToVertexOperatorId,
                                            vertexToEdgeOperatorId,
                                            edgeToEdgeOperatorId,
                                            vertexDoFDstId,
                                            vertexDoFRhsId,
                                            edgeDoFDstId,
                                            edgeDoFRhsId,
                                            backwards );
#endif
      WALBERLA_UNUSED( vertexToVertexOperatorMapId );
   }
   else
   {
      smoothSOR3DUpdateVertexDoFs( level,
                                   storage,
                                   edge,
                                   relax,
                                   vertexToVertexOperatorId,
                                   edgeToVertexOperatorId,
                                   vertexToEdgeOperatorId,
                                   edgeToEdgeOperatorId,
                                   vertexDoFDstId,
                                   vertexDoFRhsId,
                                   edgeDoFDstId,
                                   edgeDoFRhsId,
                                   backwards );
   }
   WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "VertexDoFs" ); }

   WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "EdgeDoFs" ); }
   smoothSOR3DUpdateEdgeDoFs( level,
                              storage,
                              edge,
                              relax,
                              vertexToVertexOperatorId,
                              edgeToVertexOperatorId,
                              vertexToEdgeOperatorId,
                              edgeToEdgeOperatorId,
                              vertexDoFDstId,
                              vertexDoFRhsId,
                              edgeDoFDstId,
                              edgeDoFRhsId,
                              backwards );
   WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "EdgeDoFs" ); }
}

void smoothJacobi( const uint_t&                                            level,
                   const Edge&                                              edge,
                   const real_t&                                            relax,
                   const PrimitiveDataID< StencilMemory < real_t >, Edge >& vertexToVertexStencilID,
                   const PrimitiveDataID< StencilMemory < real_t >, Edge >& edgeToVertexStencilID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcVertexDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstVertexDoFID,
                   const PrimitiveDataID< StencilMemory < real_t >, Edge >& vertexToEdgeStencilID,
                   const PrimitiveDataID< StencilMemory < real_t >, Edge >& edgeToEdgeStencilID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcEdgeDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstEdgeDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsVertexDoFID,
                   const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsEdgeDoFID )
{
   real_t* vertexToVertexStencil = edge.getData( vertexToVertexStencilID )->getPointer( level );
   real_t* edgeToVertexStencil   = edge.getData( edgeToVertexStencilID )->getPointer( level );
   real_t* srcVertexDoF          = edge.getData( srcVertexDoFID )->getPointer( level );
   real_t* dstVertexDoF          = edge.getData( dstVertexDoFID )->getPointer( level );
   real_t* vertexToEdgeStencil   = edge.getData( vertexToEdgeStencilID )->getPointer( level );
   real_t* edgeToEdgeStencil     = edge.getData( edgeToEdgeStencilID )->getPointer( level );
   real_t* srcEdgeDoF            = edge.getData( srcEdgeDoFID )->getPointer( level );
   real_t* dstEdgeDoF            = edge.getData( dstEdgeDoFID )->getPointer( level );
   real_t* rhsVertexDoF          = edge.getData( rhsVertexDoFID )->getPointer( level );
   real_t* rhsEdgeDoF            = edge.getData( rhsEdgeDoFID )->getPointer( level );

   real_t tmpVertex = 0, tmpEdgeHO = 0;

   const real_t invVertexCenter = real_c( 1.0 ) / vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
   const real_t invEdgeXCenter  = real_c( 1.0 ) / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )];

   for( const auto& it : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      ////////// VERTEX //////////
      if( it.x() != 0 )
      {
         tmpVertex = rhsVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )];
         /// on edge vertex dof
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
         {
            tmpVertex -= srcVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// on edge edge dof
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnEdgeFromVertex )
         {
            tmpVertex -= srcEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         /// south face vertex dof
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
         {
            tmpVertex -= srcVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
         }
         /// south face edge
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnSouthFaceFromVertex )
         {
            tmpVertex -= srcEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                         edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
         }
         if( edge.getNumNeighborFaces() == 2 )
         {
            /// north face vertex dof
            for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
            {
               tmpVertex -= srcVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), dir )] *
                            vertexToVertexStencil[vertexdof::stencilIndexFromVertex( dir )];
            }
            /// north face edge
            for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnNorthFaceFromVertex )
            {
               tmpVertex -= srcEdgeDoF[edgedof::macroedge::indexFromVertex( level, it.x(), dir )] *
                            edgeToVertexStencil[edgedof::stencilIndexFromVertex( dir )];
            }
         }
         dstVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )] =
             (real_c( 1.0 ) - relax) * srcVertexDoF[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )] +
             relax * invVertexCenter * tmpVertex;
      }
      ////////// HORIZONTAL EDGE //////////
      tmpEdgeHO = rhsEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )];
      /// on edge
      for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on south face
      for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF )
      {
         tmpEdgeHO -= srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
      }

      for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge )
      {
         tmpEdgeHO -= srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                      edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
      }
      /// on north face
      if( edge.getNumNeighborFaces() == 2 )
      {
         for( const auto& dir : hyteg::vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF )
         {
            tmpEdgeHO -= srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( dir )];
         }
         for( const auto& dir : hyteg::edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge )
         {
            tmpEdgeHO -= srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), dir )] *
                         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( dir )];
         }
      }
      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )] =
          (real_c( 1.0 ) - relax) * srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )] +
          relax * invEdgeXCenter * tmpEdgeHO;
   }
}

} // namespace macroedge
} // namespace P2
} // namespace hyteg
