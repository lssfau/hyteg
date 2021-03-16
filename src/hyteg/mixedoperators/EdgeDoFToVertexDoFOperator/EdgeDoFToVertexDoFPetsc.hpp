/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {
namespace EdgeDoFToVertexDoF {

using walberla::real_t;
using walberla::uint_t;

#ifdef HYTEG_BUILD_WITH_PETSC

inline void saveVertexOperator( const uint_t&                                                level,
                                const Vertex&                                                vertex,
                                const PrimitiveDataID< StencilMemory< real_t >, Vertex >&    operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex >& srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex >& dstId,
                                const std::shared_ptr< SparseMatrixProxy >&                  mat )
{
   auto opr_data = vertex.getData( operatorId )->getPointer( level );
   auto src      = vertex.getData( srcId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );

   WALBERLA_ASSERT_LESS_EQUAL( vertex.getNumNeighborEdges() + vertex.getNumNeighborFaces(),
                               vertex.getData( srcId )->getSize( level ),
                               "Stencil memory size is smaller than it should be." );

   for ( uint_t i = 0; i < vertex.getNumNeighborEdges() + vertex.getNumNeighborFaces(); i++ )
   {
      mat->addValue( uint_c( dst[0] ), uint_c( src[i] ), opr_data[i] );
   }
}

inline void saveVertexOperator3D( const uint_t&                                                                level,
                                  const Vertex&                                                                vertex,
                                  const PrimitiveStorage&                                                      storage,
                                  const PrimitiveDataID< LevelWiseMemory< MacroVertexStencilMap_T >, Vertex >& operatorId,
                                  const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex >&                 srcId,
                                  const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex >&                 dstId,
                                  const std::shared_ptr< SparseMatrixProxy >&                                  mat )
{
   auto opr_data = vertex.getData( operatorId )->getData( level );
   auto src      = vertex.getData( srcId )->getPointer( level );
   auto dst      = vertex.getData( dstId )->getPointer( level );

   const auto centerIndexOnVertex = indexing::Index( 0, 0, 0 );

   for ( uint_t neighborCellID = 0; neighborCellID < vertex.getNumNeighborCells(); neighborCellID++ )
   {
      const Cell& neighborCell      = *( storage.getCell( vertex.neighborCells().at( neighborCellID ) ) );
      auto        cellLocalVertexID = neighborCell.getLocalVertexID( vertex.getID() );

      const auto basisInCell = algorithms::getMissingIntegersAscending< 1, 4 >( {cellLocalVertexID} );

      const auto centerIndexInCell = indexing::basisConversion(
          centerIndexOnVertex, basisInCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );

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
               WALBERLA_ASSERT_GREATER_EQUAL( level, 1 );
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

            const auto dstInt = dst[0];
            const auto srcInt = src[leafArrayIndexOnVertex];
            mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
         }
      }
   }
}

inline void saveEdgeOperator( const uint_t&                                              Level,
                              const Edge&                                                edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge >&    operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge >& srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge >& dstId,
                              const std::shared_ptr< SparseMatrixProxy >&                mat )
{
   const real_t*   opr_data = edge.getData( operatorId )->getPointer( Level );
   const PetscInt* src      = edge.getData( srcId )->getPointer( Level );
   const PetscInt* dst      = edge.getData( dstId )->getPointer( Level );

   PetscInt srcInt;
   PetscInt dstInt;

   for ( const auto it : vertexdof::macroedge::Iterator( Level, 1 ) )
   {
      dstInt = dst[vertexdof::macroedge::indexFromVertex( Level, it.col(), stencilDirection::VERTEX_C )];

      for ( const auto& neighbor : edgedof::macroedge::neighborsOnEdgeFromVertex )
      {
         srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
         mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[edgedof::stencilIndexFromVertex( neighbor )] );
      }

      for ( const auto& neighbor : edgedof::macroedge::neighborsOnSouthFaceFromVertex )
      {
         srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
         mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[edgedof::stencilIndexFromVertex( neighbor )] );
      }

      if ( edge.getNumNeighborFaces() == 2 )
      {
         for ( const auto& neighbor : edgedof::macroedge::neighborsOnNorthFaceFromVertex )
         {
            srcInt = src[edgedof::macroedge::indexFromVertex( Level, it.col(), neighbor )];
            mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[edgedof::stencilIndexFromVertex( neighbor )] );
         }
      }
   }
}

inline void saveEdgeOperator3D( const uint_t&                                                            level,
                                const Edge&                                                              edge,
                                const PrimitiveStorage&                                                  storage,
                                const PrimitiveDataID< LevelWiseMemory< MacroEdgeStencilMap_T >, Edge >& operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge >&               srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Edge >&               dstId,
                                const std::shared_ptr< SparseMatrixProxy >&                              mat )
{
   auto opr_data = edge.getData( operatorId )->getData( level );
   auto src      = edge.getData( srcId )->getPointer( level );
   auto dst      = edge.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexOnEdge : hyteg::vertexdof::macroedge::Iterator( level, 1 ) )
   {
      const auto dstInt = dst[vertexdof::macroedge::index( level, centerIndexOnEdge.x() )];

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
              neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

         const auto centerIndexInCell = indexing::basisConversion(
             centerIndexOnEdge, basisInCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );

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

inline void saveFaceOperator( const uint_t&                                              Level,
                              const Face&                                                face,
                              const PrimitiveDataID< StencilMemory< real_t >, Face >&    operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face >& srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face >& dstId,
                              const std::shared_ptr< SparseMatrixProxy >&                mat )
{
   const real_t*   opr_data = face.getData( operatorId )->getPointer( Level );
   const PetscInt* src      = face.getData( srcId )->getPointer( Level );
   const PetscInt* dst      = face.getData( dstId )->getPointer( Level );

   PetscInt srcInt;
   PetscInt dstInt;

   for ( const auto& it : vertexdof::macroface::Iterator( Level, 1 ) )
   {
      dstInt = dst[vertexdof::macroface::indexFromVertex( Level, it.col(), it.row(), stencilDirection::VERTEX_C )];

      for ( const auto& neighbor : edgedof::macroface::neighborsFromVertex )
      {
         srcInt = src[edgedof::macroface::indexFromVertex( Level, it.col(), it.row(), neighbor )];
         mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[edgedof::stencilIndexFromVertex( neighbor )] );
      }
   }
}

inline void saveFaceOperator3D( const uint_t&                                                            level,
                                const Face&                                                              face,
                                const PrimitiveStorage&                                                  storage,
                                const PrimitiveDataID< LevelWiseMemory< MacroFaceStencilMap_T >, Face >& operatorId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face >&               srcId,
                                const PrimitiveDataID< FunctionMemory< PetscInt >, Face >&               dstId,
                                const std::shared_ptr< SparseMatrixProxy >&                              mat )
{
   auto opr_data = face.getData( operatorId )->getData( level );
   auto src      = face.getData( srcId )->getPointer( level );
   auto dst      = face.getData( dstId )->getPointer( level );

   for ( const auto& centerIndexInFace : hyteg::vertexdof::macroface::Iterator( level, 1 ) )
   {
      const auto dstIdx = vertexdof::macroface::index( level, centerIndexInFace.x(), centerIndexInFace.y() );
      const auto dstInt = dst[dstIdx];

      for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
      {
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

               WALBERLA_ASSERT_LESS( leafArrayIndexInFace, face.getData( srcId )->getSize( level ) );
               const auto srcInt = src[leafArrayIndexInFace];
               mat->addValue( uint_c( dstInt ), uint_c( srcInt ), stencilWeight );
            }
         }
      }
   }
}

inline void
    saveCellOperator( const uint_t&                                                                                Level,
                      const Cell&                                                                                  cell,
                      const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >& operatorId,
                      const PrimitiveDataID< FunctionMemory< PetscInt >, Cell >&                                   srcId,
                      const PrimitiveDataID< FunctionMemory< PetscInt >, Cell >&                                   dstId,
                      const std::shared_ptr< SparseMatrixProxy >&                                                  mat )
{
   auto      opr_data = cell.getData( operatorId )->getData( Level );
   PetscInt* src      = cell.getData( srcId )->getPointer( Level );
   PetscInt* dst      = cell.getData( dstId )->getPointer( Level );

   for ( const auto& it : vertexdof::macrocell::Iterator( Level, 1 ) )
   {
      const auto dstArrayIdx = vertexdof::macrocell::index( Level, it.x(), it.y(), it.z() );
      const auto dstInt      = dst[dstArrayIdx];

      for ( const auto& orientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromVertexDoFInMacroCell( orientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            const auto srcIdx      = it + neighbor;
            const auto srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), orientation );
            const auto srcInt      = src[srcArrayIdx];
            mat->addValue( uint_c( dstInt ), uint_c( srcInt ), opr_data[orientation][neighbor] );
         }
      }
   }
}

} // namespace EdgeDoFToVertexDoF

namespace petsc {

template < class OperatorType >
inline void createMatrix( const OperatorType&                         opr,
                          const EdgeDoFFunction< PetscInt >&          src,
                          const P1Function< PetscInt >&               dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          size_t                                      level,
                          DoFType                                     flag )
{
   const auto storage = src.getStorage();

   for ( auto& it : opr.getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage->hasGlobalCells() )
         {
            EdgeDoFToVertexDoF::saveVertexOperator3D(
                level, vertex, *storage, opr.getVertexStencil3DID(), src.getVertexDataID(), dst.getVertexDataID(), mat );
         }
         else
         {
            EdgeDoFToVertexDoF::saveVertexOperator(
                level, vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat );
         }
      }
   }

   if ( level >= 1 )
   {
      for ( auto& it : opr.getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            if ( storage->hasGlobalCells() )
            {
               EdgeDoFToVertexDoF::saveEdgeOperator3D(
                   level, edge, *storage, opr.getEdgeStencil3DID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
            }
            else
            {
               EdgeDoFToVertexDoF::saveEdgeOperator(
                   level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
            }
         }
      }
   }

   if ( level >= 2 )
   {
      for ( auto& it : opr.getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage->hasGlobalCells() )
            {
               EdgeDoFToVertexDoF::saveFaceOperator3D(
                   level, face, *storage, opr.getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
            else
            {
               EdgeDoFToVertexDoF::saveFaceOperator(
                   level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
         }
      }

      for ( auto& it : opr.getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            EdgeDoFToVertexDoF::saveCellOperator(
                level, cell, opr.getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
         }
      }
   }
}

#endif

} // namespace EdgeDoFToVertexDoF

} // namespace hyteg
