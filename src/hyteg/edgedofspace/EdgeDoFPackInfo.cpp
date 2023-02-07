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
#include "EdgeDoFPackInfo.hpp"

#include "hyteg/Algorithms.hpp"
#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/generatedKernels/communicate_buffered_pack_edgedof_face_to_cell.hpp"
#include "hyteg/edgedofspace/generatedKernels/communicate_buffered_unpack_edgedof_face_to_cell.hpp"
#include "hyteg/edgedofspace/generatedKernels/communicate_directly_edgedof_cell_to_face_part_1.hpp"
#include "hyteg/edgedofspace/generatedKernels/communicate_directly_edgedof_cell_to_face_part_2.hpp"
#include "hyteg/edgedofspace/generatedKernels/communicate_directly_edgedof_face_to_cell.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {

template < typename ValueType >
EdgeDoFPackInfo< ValueType >::EdgeDoFPackInfo( uint_t                                                 level,
                                               PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                               PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                               PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                               PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                                               std::weak_ptr< PrimitiveStorage >                      storage )
: communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
{}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   ValueType* edgeData       = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     vertexIdOnEdge = sender->vertex_index( receiver );
   if( vertexIdOnEdge == 0 )
   {
      buffer << edgeData[edgedof::macroedge::index( level_, 0 )];
      for( const PrimitiveID& faceID : sender->neighborFaces() )
      {
         buffer << edgeData[edgedof::macroedge::indexOnNeighborFace(
             level_, 0, sender->face_index( faceID ), edgedof::EdgeDoFOrientation::XY )];
      }
   } else if( vertexIdOnEdge == 1 )
   {
      uint_t length = levelinfo::num_microedges_per_edge( level_ );
      buffer << edgeData[edgedof::macroedge::index( level_, idx_t( length - 1 ) )];
      for( const PrimitiveID& faceID : sender->neighborFaces() )
      {
         buffer << edgeData[edgedof::macroedge::indexOnNeighborFace(
             level_, idx_t( length - 1 ), sender->face_index( faceID ), edgedof::EdgeDoFOrientation::Y )];
      }
   } else
   {
      WALBERLA_ABORT( "vertex is not part of edge" )
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   buffer >> vertexData[receiver->edge_index( sender )];
   for( const PrimitiveID& faceID : storage_.lock()->getEdge( sender )->neighborFaces() )
   {
      buffer >> vertexData[receiver->getNumNeighborEdges() + receiver->face_index( faceID )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   ValueType* edgeData       = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     vertexIdOnEdge = sender->vertex_index( receiver->getID() );
   ValueType* vertexData     = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   uint_t     edgeIdOnVertex = receiver->edge_index( sender->getID() );

   if( vertexIdOnEdge == 0 )
   {
      vertexData[edgeIdOnVertex] = edgeData[edgedof::macroedge::index( level_, 0 )];
      for( const PrimitiveID& faceID : sender->neighborFaces() )
      {
         vertexData[receiver->getNumNeighborEdges() + receiver->face_index( faceID )] =
             edgeData[edgedof::macroedge::indexOnNeighborFace(
                 level_, 0, sender->face_index( faceID ), edgedof::EdgeDoFOrientation::XY )];
      }
   } else if( vertexIdOnEdge == 1 )
   {
      uint_t edgeLength          = levelinfo::num_microedges_per_edge( level_ );
      vertexData[edgeIdOnVertex] = edgeData[edgedof::macroedge::index( level_, idx_t( edgeLength - 1 ) )];
      for( const PrimitiveID& faceID : sender->neighborFaces() )
      {
         vertexData[receiver->getNumNeighborEdges() + receiver->face_index( faceID )] =
             edgeData[edgedof::macroedge::indexOnNeighborFace(
                 level_, idx_t( edgeLength - 1 ), sender->face_index( faceID ), edgedof::EdgeDoFOrientation::Y )];
      }
   } else
   {
      WALBERLA_ABORT( "vertex is not part of edge" )
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_UNUSED( receiver );
   ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
   for( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); ++i )
   {
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, idx_t( i ), stencilDirection::EDGE_HO_C )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   using edgedof::macroface::indexFromHorizontalEdge;
   using hyteg::edgedof::macroface::BoundaryIterator;
   ValueType*                    faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = receiver->edge_index( sender );
   indexing::FaceBoundaryDirection faceDir =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );
   for( const auto& it : BoundaryIterator( level_, faceDir, 0 ) )
   {
      if( edgeIndexOnFace == 0 )
      {
         buffer >>
             faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 2 )
      {
         buffer >>
             faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_DI_N )];
      } else if( edgeIndexOnFace == 1 )
      {
         buffer >>
             faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_VE_NW )];
      } else
      {
         WALBERLA_ABORT( "Wrong edgeIndexOnFace" )
      }
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Edge to Face" );
   using edgedof::macroface::indexFromHorizontalEdge;
   using hyteg::edgedof::macroface::BoundaryIterator;
   ValueType*                    faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   ValueType*                    edgeData        = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = receiver->edge_index( sender->getID() );
   indexing::FaceBoundaryDirection faceDir =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );
   uint_t indexOnEdge = 0;
   for( const auto& it : BoundaryIterator( level_, faceDir, 0 ) )
   {
      if( edgeIndexOnFace == 0 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_HO_C )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, idx_t( indexOnEdge ), stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 2 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_DI_N )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, idx_t( indexOnEdge ), stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 1 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_VE_NW )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, idx_t( indexOnEdge ), stencilDirection::EDGE_HO_C )];
      } else
      {
         WALBERLA_ABORT( "Wrong edgeIndexOnFace" )
      }
      ++indexOnEdge;
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Edge to Face" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Edge (pack)" );
   using hyteg::edgedof::macroface::BoundaryIterator;
   ValueType*                    faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = sender->edge_index( receiver );
   indexing::FaceBoundaryDirection faceBorderDir =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, sender->getEdgeOrientation()[edgeIndexOnFace] );
   stencilDirection faceDirOne;
   stencilDirection faceDirTwo;
   stencilDirection faceDirThree;
   if( edgeIndexOnFace == 0 )
   {
      faceDirOne = stencilDirection::EDGE_HO_C;
      if( sender->getEdgeOrientation()[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_DI_N;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_DI_N;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      }
   } else if( edgeIndexOnFace == 2 )
   {
      faceDirOne = stencilDirection::EDGE_DI_N;
      if( sender->getEdgeOrientation()[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_HO_C;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_HO_C;
      }
   } else if( edgeIndexOnFace == 1 )
   {
      faceDirOne = stencilDirection::EDGE_VE_NW;
      if( sender->getEdgeOrientation()[edgeIndexOnFace] == -1 )
      {
         faceDirTwo   = stencilDirection::EDGE_DI_N;
         faceDirThree = stencilDirection::EDGE_HO_C;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_HO_C;
         faceDirThree = stencilDirection::EDGE_DI_N;
      }
   } else
   {
      WALBERLA_ABORT( "Wrong edgeIndexOnFace" )
   }
   for( const auto& it : BoundaryIterator( level_, faceBorderDir, 1 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirOne )];
   }
   for( const auto& it : BoundaryIterator( level_, faceBorderDir, 0 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirTwo )];
   }
   for( const auto& it : BoundaryIterator( level_, faceBorderDir, 0 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirThree )];
   }

   //// DoFs on neighboring cells /////
   for ( const auto& neighborCellID : sender->neighborCells() )
   {
      const Cell&  neighborCell    = *( storage_.lock()->getCell( neighborCellID ) );
      const auto   cellLocalEdgeID = neighborCell.getLocalEdgeID( receiver );
      const uint_t faceLocalCellID = sender->cell_index( neighborCellID );

      const uint_t cellLocalVertexID0OfFace = neighborCell.getLocalVertexID( sender->getVertexID0() );
      const uint_t cellLocalVertexID1OfFace = neighborCell.getLocalVertexID( sender->getVertexID1() );
      const uint_t cellLocalVertexID2OfFace = neighborCell.getLocalVertexID( sender->getVertexID2() );
      const uint_t cellLocalRemainingVertex =
          6 - ( cellLocalVertexID0OfFace + cellLocalVertexID1OfFace + cellLocalVertexID2OfFace );

      const std::array< uint_t, 4 > faceBasisInCell = {
          cellLocalVertexID0OfFace, cellLocalVertexID1OfFace, cellLocalVertexID2OfFace, cellLocalRemainingVertex};

      const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
          {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
           neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

      for ( const auto edgeOrientationOnReferenceEdge : edgedof::allEdgeDoFOrientations )
      {
         if ( level_ == 0 && edgeOrientationOnReferenceEdge != edgedof::EdgeDoFOrientation::YZ )
         {
            continue;
         }

         for ( const auto& indexOnEdge : hyteg::edgedof::macroedge::Iterator( level_, 0 ) )
         {
            auto indexOnEdgeCopy = indexOnEdge;

            switch ( edgeOrientationOnReferenceEdge )
            {
            case edgedof::EdgeDoFOrientation::X:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 2 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               indexOnEdgeCopy.y()++;
               break;
            case edgedof::EdgeDoFOrientation::Y:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               break;
            case edgedof::EdgeDoFOrientation::Z:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.y()++;
               break;
            case edgedof::EdgeDoFOrientation::XY:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               break;
            case edgedof::EdgeDoFOrientation::XZ:
               indexOnEdgeCopy.y()++;
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               break;
            case edgedof::EdgeDoFOrientation::YZ:
               break;
            case edgedof::EdgeDoFOrientation::XYZ:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               break;
            default:
               WALBERLA_ABORT( "wrong direction" );
            }

            /// the size has to be adjusted for the XYZ edge
            uint_t cellWidth;
            if ( edgeOrientationOnReferenceEdge == edgedof::EdgeDoFOrientation::XYZ )
            {
               cellWidth = levelinfo::num_microedges_per_edge( level_ ) - 1;
            }
            else
            {
               cellWidth = levelinfo::num_microedges_per_edge( level_ );
            }
            const auto indexInCell           = indexing::basisConversion( indexOnEdgeCopy, basisInCell, {0, 1, 2, 3}, cellWidth );
            const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                edgeOrientationOnReferenceEdge, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

            const auto indexOnFace           = indexing::basisConversion( indexInCell, {0, 1, 2, 3}, faceBasisInCell, cellWidth );
            const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationCellToFace(
                cellCenterOrientation, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ), faceBasisInCell.at( 2 ) );

            uint_t idxOnFace =
                edgedof::macroface::index( level_, indexOnFace.x(), indexOnFace.y(), edgeOrientationOnFace, faceLocalCellID );
            buffer << faceData[idxOnFace];
         }
      }
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Edge (pack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Edge (unpack)" );
   ValueType* edgeData        = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     edgeLocalFaceID = receiver->face_index( sender );
   /////////// DoFs on Face ///////////
   for( const auto edgeOrienation : edgedof::faceLocalEdgeDoFOrientations )
   {
      for( const auto& indexOnEdge : hyteg::edgedof::macroedge::Iterator( level_, 0 ) )
      {
         if( edgeOrienation == edgedof::EdgeDoFOrientation::X &&
             indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
         {
            continue;
         }
         uint_t idxOnEdge = edgedof::macroedge::indexOnNeighborFace( level_, indexOnEdge.x(), edgeLocalFaceID, edgeOrienation );
         buffer >> edgeData[idxOnEdge];
      }
   }

   //// DoFs on neighboring cells /////
   for ( const auto& neighborCellID : storage_.lock()->getFace( sender )->neighborCells() )
   {
      const uint_t edgeLocalCellID = receiver->cell_index( neighborCellID );

      for ( const auto edgeOrientation : edgedof::allEdgeDoFOrientations )
      {
         for ( const auto& indexOnEdge : hyteg::edgedof::macroedge::Iterator( level_, 0 ) )
         {
            if ( level_ == 0 && edgeOrientation != edgedof::EdgeDoFOrientation::YZ )
            {
               continue;
            }

            switch ( edgeOrientation )
            {
            case edgedof::EdgeDoFOrientation::X:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 2 ) )
                  continue;
               break;
            case edgedof::EdgeDoFOrientation::YZ:
               break;
            default:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               break;
            }

            uint_t idxOnEdge = edgedof::macroedge::indexOnNeighborCell(
                level_, indexOnEdge.x(), edgeLocalCellID, receiver->getNumNeighborFaces(), edgeOrientation );
            buffer >> edgeData[idxOnEdge];
         }
      }
   }

   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Edge (unpack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Edge" );
   ValueType* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   uint_t edgeIndexCounter;
   uint_t edgeLocalFaceID          = receiver->face_index( sender->getID() );
   uint_t faceLocalEdgeID          = sender->edge_index( receiver->getID() );
   uint_t faceLocalVertexIDOfEdge0 = sender->vertex_index( receiver->getVertexID0() );
   uint_t faceLocalVertexIDOfEdge1 = sender->vertex_index( receiver->getVertexID1() );

   indexing::FaceBoundaryDirection faceBorderDir =
       indexing::getFaceBoundaryDirection( faceLocalEdgeID, sender->getEdgeOrientation()[faceLocalEdgeID] );

   /////////// DoFs on Face ///////////
   for( const auto edgeOriOnReferenceEdge : {hyteg::edgedof::EdgeDoFOrientation::Y, hyteg::edgedof::EdgeDoFOrientation::XY} )
   {
      edgeIndexCounter = 0;

      const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationEdgeToFace(
          edgeOriOnReferenceEdge, faceLocalVertexIDOfEdge0, faceLocalVertexIDOfEdge1 );

      for( const auto& it : edgedof::macroface::BoundaryIterator( level_, faceBorderDir, 0 ) )
      {
         uint_t idxOnEdge =
             edgedof::macroedge::indexOnNeighborFace( level_, idx_t( edgeIndexCounter ), edgeLocalFaceID, edgeOriOnReferenceEdge );
         uint_t idxOnFace    = edgedof::macroface::index( level_, it.col(), it.row(), edgeOrientationOnFace );
         edgeData[idxOnEdge] = faceData[idxOnFace];
         ++edgeIndexCounter;
      }
   }
   for( const auto edgeOriOnReferenceEdge : {hyteg::edgedof::EdgeDoFOrientation::X} )
   {
      edgeIndexCounter                 = 0;
      const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationEdgeToFace(
          edgeOriOnReferenceEdge, faceLocalVertexIDOfEdge0, faceLocalVertexIDOfEdge1 );

      for( const auto& it : edgedof::macroface::BoundaryIterator( level_, faceBorderDir, 1 ) )
      {
         uint_t idxOnEdge =
             edgedof::macroedge::indexOnNeighborFace( level_, idx_t( edgeIndexCounter ), edgeLocalFaceID, edgeOriOnReferenceEdge );
         edgeData[idxOnEdge] = faceData[edgedof::macroface::index( level_, it.col(), it.row(), edgeOrientationOnFace )];
         ++edgeIndexCounter;
      }
   }
   ////////////////////////////////////

   //// DoFs on neighboring cells /////
   for ( const auto& neighborCellID : sender->neighborCells() )
   {
      const Cell&  neighborCell    = *( storage_.lock()->getCell( neighborCellID ) );
      const auto   cellLocalEdgeID = neighborCell.getLocalEdgeID( receiver->getID() );
      const uint_t edgeLocalCellID = receiver->cell_index( neighborCellID );
      const uint_t faceLocalCellID = sender->cell_index( neighborCellID );

      const uint_t cellLocalVertexID0OfFace = neighborCell.getLocalVertexID( sender->getVertexID0() );
      const uint_t cellLocalVertexID1OfFace = neighborCell.getLocalVertexID( sender->getVertexID1() );
      const uint_t cellLocalVertexID2OfFace = neighborCell.getLocalVertexID( sender->getVertexID2() );
      const uint_t cellLocalRemainingVertex =
          6 - ( cellLocalVertexID0OfFace + cellLocalVertexID1OfFace + cellLocalVertexID2OfFace );

      const std::array< uint_t, 4 > faceBasisInCell = {
          cellLocalVertexID0OfFace, cellLocalVertexID1OfFace, cellLocalVertexID2OfFace, cellLocalRemainingVertex};

      const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
          {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
           neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

      for ( const auto edgeOrientationOnReferenceEdge : edgedof::allEdgeDoFOrientations )
      {
         if ( level_ == 0 && edgeOrientationOnReferenceEdge != edgedof::EdgeDoFOrientation::YZ )
         {
           continue;
         }

         for ( const auto& indexOnEdge : hyteg::edgedof::macroedge::Iterator( level_, 0 ) )
         {
            auto indexOnEdgeCopy = indexOnEdge;

            switch ( edgeOrientationOnReferenceEdge )
            {
            case edgedof::EdgeDoFOrientation::X:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 2 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               indexOnEdgeCopy.y()++;
               break;
            case edgedof::EdgeDoFOrientation::Y:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               break;
            case edgedof::EdgeDoFOrientation::Z:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.y()++;
               break;
            case edgedof::EdgeDoFOrientation::XY:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               indexOnEdgeCopy.z()++;
               break;
            case edgedof::EdgeDoFOrientation::XZ:
               indexOnEdgeCopy.y()++;
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               break;
            case edgedof::EdgeDoFOrientation::YZ:
               break;
            case edgedof::EdgeDoFOrientation::XYZ:
               if ( indexOnEdge.x() >= idx_t( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
                  continue;
               break;
            default:
               WALBERLA_ABORT( "wrong direction" );
            }

            /// the size has to be adjusted for the XYZ edge
            uint_t cellWidth;
            if ( edgeOrientationOnReferenceEdge == edgedof::EdgeDoFOrientation::XYZ )
            {
               cellWidth = levelinfo::num_microedges_per_edge( level_ ) - 1;
            }
            else
            {
               cellWidth = levelinfo::num_microedges_per_edge( level_ );
            }
            const auto indexInCell           = indexing::basisConversion( indexOnEdgeCopy, basisInCell, {0, 1, 2, 3}, cellWidth );
            const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                edgeOrientationOnReferenceEdge, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

            const auto indexOnFace           = indexing::basisConversion( indexInCell, {0, 1, 2, 3}, faceBasisInCell, cellWidth );
            const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationCellToFace(
                cellCenterOrientation, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ), faceBasisInCell.at( 2 ) );

            uint_t idxOnEdge = edgedof::macroedge::indexOnNeighborCell(
                level_, indexOnEdge.x(), edgeLocalCellID, receiver->getNumNeighborFaces(), edgeOrientationOnReferenceEdge );
            uint_t idxOnFace =
                edgedof::macroface::index( level_, indexOnFace.x(), indexOnFace.y(), edgeOrientationOnFace, faceLocalCellID );
            edgeData[idxOnEdge] = faceData[idxOnFace];
         }
      }
   }

   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Edge" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell (pack)" );
   WALBERLA_UNUSED( receiver );
   const ValueType* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
   for( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell (pack)" );
}

template <>
void EdgeDoFPackInfo< real_t >::packFaceForCell( const Face*                sender,
                                                 const PrimitiveID&         receiver,
                                                 walberla::mpi::SendBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell (pack)" );
   if ( globalDefines::useGeneratedKernels && level_ >= 1 )
   {
#ifdef HYTEG_USE_GENERATED_KERNELS
      auto         cell             = storage_.lock()->getCell( receiver );
      const uint_t localFaceID      = cell->getLocalFaceID( sender->getID() );
      const uint_t iterationVertex0 = cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
      const uint_t iterationVertex1 = cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
      const uint_t iterationVertex2 = cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdx;
      for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
         firstIdx[e] = edgedof::macroface::index( level_, 0, 0, e );
      const uint_t  requiredBufferElements = levelinfo::num_microedges_per_face( level_ );
      auto          buffer_ptr             = (real_t*) ( buffer.forward( requiredBufferElements * sizeof( real_t ) ) );
      const real_t* faceData               = sender->getData( dataIDFace_ )->getPointer( level_ );
      edgedof::comm::generated::communicate_buffered_pack_edgedof_face_to_cell( &faceData[firstIdx[eo::X]],
                                                                                &faceData[firstIdx[eo::XY]],
                                                                                &faceData[firstIdx[eo::Y]],
                                                                                buffer_ptr,
                                                                                static_cast< int32_t >( level_ ),
                                                                                static_cast< int >( iterationVertex0 ),
                                                                                static_cast< int >( iterationVertex1 ),
                                                                                static_cast< int >( iterationVertex2 ),
                                                                                0 );
#endif
   }
   else
   {
      WALBERLA_UNUSED( receiver );
      const real_t* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
      for ( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
      {
         buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
         buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
         buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];
      }
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell (pack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell (unpack)" );
   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = receiver->getLocalFaceID( sender );
   const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   for( const auto& cellIterator :
        edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      buffer >> cellData[edgedof::macrocell::index(
                    level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationX )];
      buffer >> cellData[edgedof::macrocell::index(
                    level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationY )];
      buffer >> cellData[edgedof::macrocell::index(
                    level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationXY )];
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell (unpack)" );
}

template <>
void EdgeDoFPackInfo< real_t >::unpackCellFromFace( Cell*                      receiver,
                                                    const PrimitiveID&         sender,
                                                    walberla::mpi::RecvBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell (unpack)" );
   real_t* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = receiver->getLocalFaceID( sender );
   const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   if ( globalDefines::useGeneratedKernels && level_ >= 1 )
   {
#ifdef HYTEG_USE_GENERATED_KERNELS
      const uint_t requiredBufferElements = levelinfo::num_microedges_per_face( level_ );
      auto         buffer_ptr             = (real_t*) ( buffer.skip( requiredBufferElements * sizeof( real_t ) ) );

      std::map< edgedof::EdgeDoFOrientation, uint_t > firstIdxCell;
      for ( auto e : edgedof::allEdgeDoFOrientations )
         firstIdxCell[e] = edgedof::macrocell::index( level_, 0, 0, 0, e );

      edgedof::comm::generated::communicate_buffered_unpack_edgedof_face_to_cell(
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::X]],
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XY]],
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XZ]],
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Y]],
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::YZ]],
          &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Z]],
          buffer_ptr,
          static_cast< int32_t >( level_ ),
          static_cast< int >( iterationVertex0 ),
          static_cast< int >( iterationVertex1 ),
          static_cast< int >( iterationVertex2 ),
          0 );
#endif
   }
   else
   {
      auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
          edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
      auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
          edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
      auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
          edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

      for ( const auto& cellIterator :
            edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
      {
         buffer >> cellData[edgedof::macrocell::index(
                       level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationX )];
         buffer >> cellData[edgedof::macrocell::index(
                       level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationY )];
         buffer >> cellData[edgedof::macrocell::index(
                       level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationXY )];
      }
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell (unpack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell" );
   const ValueType* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
   ValueType*       cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = receiver->getLocalFaceID( sender->getID() );
   const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;

      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationX )] =
          faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationY )] =
          faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationXY )] =
          faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];

      cellIterator++;
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell" );
}


template <>
void EdgeDoFPackInfo< real_t >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
  this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Face to Cell" );
  const real_t* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
  real_t*       cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

  const uint_t localFaceID      = receiver->getLocalFaceID( sender->getID() );
  const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  if ( globalDefines::useGeneratedKernels && level_ >= 1 )
  {
#ifdef HYTEG_USE_GENERATED_KERNELS
    std::map< edgedof::EdgeDoFOrientation , uint_t > firstIdxFace;
    for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      firstIdxFace[e] = edgedof::macroface::index( level_, 0, 0, e );

    std::map< edgedof::EdgeDoFOrientation , uint_t > firstIdxCell;
    for ( auto e : edgedof::allEdgeDoFOrientations )
      firstIdxCell[e] = edgedof::macrocell::index( level_, 0, 0, 0, e );

    edgedof::comm::generated::communicate_directly_edgedof_face_to_cell(
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::X]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XY]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Y]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::YZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Z]],
      &faceData[firstIdxFace[edgedof::EdgeDoFOrientation::X]],
      &faceData[firstIdxFace[edgedof::EdgeDoFOrientation::XY]],
      &faceData[firstIdxFace[edgedof::EdgeDoFOrientation::Y]],
      static_cast< int32_t >( level_ ),
      static_cast< int >( iterationVertex0 ),
      static_cast< int >( iterationVertex1 ),
      static_cast< int >( iterationVertex2 )
    );
#endif
  }
  else
  {
    auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

    auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

    for ( const auto & faceIdx : edgedof::macroface::Iterator( level_ ))
    {
      auto cellIdx = *cellIterator;

      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationX )] =
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationY )] =
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationXY )] =
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];

      cellIterator++;
    }

    WALBERLA_ASSERT( cellIterator == cellIterator.end());
  }
  this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell" );
}


template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Cell to Face (pack)" );
   const ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto dstEdgeOrientationZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Z, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationYZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::YZ, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XZ, iterationVertex0, iterationVertex1, iterationVertex2 );

   for( const auto& cellIt :
        edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationZ )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationYZ )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationXZ )];
   }

   if ( level_ >= 1 )
   {
     auto cellItXYZ = edgedof::macrocell::BoundaryIteratorXYZ( level_, iterationVertex0, iterationVertex1, iterationVertex2 );
     for ( const auto & cellIt :
     edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 ))
     {
       buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationX )];
       buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationY )];
       buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationXY )];

       buffer << cellData[edgedof::macrocell::index(
       level_, cellItXYZ->x(), cellItXYZ->y(), cellItXYZ->z(), edgedof::EdgeDoFOrientation::XYZ )];
       cellItXYZ++;
     }
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Cell to Face (pack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Cell to Face (unpack)" );
   ValueType*   faceData    = receiver->getData( dataIDFace_ )->getPointer( level_ );
   const uint_t localCellID = receiver->cell_index( sender );

   typedef edgedof::EdgeDoFOrientation EO;

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

   for( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::Z, localCellID )];
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::YZ, localCellID )];
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::XZ, localCellID )];
   }

   if ( level_ >= 1 )
   {
     for ( const auto & faceIdx : indexing::FaceIterator( levelinfo::num_microedges_per_edge( level_ ) - 1 ))
     {
       buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::X, localCellID )];
       buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::Y, localCellID )];
       buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::XY, localCellID )];

       buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::XYZ, localCellID )];
     }
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Cell to Face (unpack)" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Cell to Face" );
   const ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
   ValueType*       faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver->getID() );
   const uint_t localCellID      = receiver->cell_index( sender->getID() );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto dstEdgeOrientationZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Z, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationYZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::YZ, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXZ = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XZ, iterationVertex0, iterationVertex1, iterationVertex2 );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;

      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Z, localCellID )] =
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationZ )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::YZ, localCellID )] =
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationYZ )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XZ, localCellID )] =
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationXZ )];

      cellIterator++;
   }

   if ( level_ >= 1 )
   {
     auto cellIterator2 = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );
     auto cellIterator3 = edgedof::macrocell::BoundaryIteratorXYZ( level_, iterationVertex0, iterationVertex1, iterationVertex2 );
     for ( const auto & faceIdx : indexing::FaceIterator( levelinfo::num_microedges_per_edge( level_ ) - 1 ))
     {
       auto cellIdx2 = *cellIterator2;

       faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X, localCellID )] =
       cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationX )];
       faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y, localCellID )] =
       cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationY )];
       faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY, localCellID )] =
       cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationXY )];

       cellIterator2++;

       auto cellIdx3 = *cellIterator3;

       faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XYZ, localCellID )] =
       cellData[edgedof::macrocell::index(
       level_, cellIdx3.x(), cellIdx3.y(), cellIdx3.z(), edgedof::EdgeDoFOrientation::XYZ )];

       cellIterator3++;
     }
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Cell to Face" );
}


template <>
void EdgeDoFPackInfo< real_t >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
  this->storage_.lock()->getTimingTree()->start( "EdgeDoF - Cell to Face" );
  const real_t * cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
  real_t *       faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

  const uint_t localFaceID      = sender->getLocalFaceID( receiver->getID() );
  const uint_t localCellID      = receiver->cell_index( sender->getID() );
  const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  if ( level_ >= 1 && globalDefines::useGeneratedKernels )
  {
#ifdef HYTEG_USE_GENERATED_KERNELS
    std::map< edgedof::EdgeDoFOrientation , uint_t > firstIdxCell;
    for ( auto e : edgedof::allEdgeDoFOrientations )
      firstIdxCell[e] = edgedof::macrocell::index( level_, 0, 0, 0, e );

    std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, uint_t > > offset_gl_orientation;
    for ( uint_t gl = 0; gl < 2; gl++ )
    {
      for ( const auto& eo : edgedof::allEdgeDoFOrientations )
      {
        offset_gl_orientation[gl][eo] = edgedof::macroface::index( level_, 0, 0, eo, gl );
      }
    }

    edgedof::comm::generated::communicate_directly_edgedof_cell_to_face_part_1(
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::X]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XY]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XYZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Y]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::YZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Z]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::X]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::XY]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::XYZ]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::Y]],
      static_cast< int32_t >( level_ ),
      static_cast< int >( iterationVertex0 ),
      static_cast< int >( iterationVertex1 ),
      static_cast< int >( iterationVertex2 )
    );

    edgedof::comm::generated::communicate_directly_edgedof_cell_to_face_part_2(
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::X]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XY]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::XZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Y]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::YZ]],
      &cellData[firstIdxCell[edgedof::EdgeDoFOrientation::Z]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::XZ]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::YZ]],
      &faceData[offset_gl_orientation[localCellID][edgedof::EdgeDoFOrientation::Z]],
      static_cast< int32_t >( level_ ),
      static_cast< int >( iterationVertex0 ),
      static_cast< int >( iterationVertex1 ),
      static_cast< int >( iterationVertex2 )
    );
#endif
  }
  else
  {
    auto dstEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

    auto dstEdgeOrientationZ = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::Z, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationYZ = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::YZ, iterationVertex0, iterationVertex1, iterationVertex2 );
    auto dstEdgeOrientationXZ = edgedof::convertEdgeDoFOrientationFaceToCell(
    edgedof::EdgeDoFOrientation::XZ, iterationVertex0, iterationVertex1, iterationVertex2 );

    WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
    WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID()));

    auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

    for ( const auto & faceIdx : edgedof::macroface::Iterator( level_ ))
    {
      auto cellIdx = *cellIterator;

      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Z, localCellID )] =
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationZ )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::YZ, localCellID )] =
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationYZ )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XZ, localCellID )] =
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationXZ )];

      cellIterator++;
    }

    if ( level_ >= 1 )
    {
      auto cellIterator2 = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );
      auto cellIterator3 = edgedof::macrocell::BoundaryIteratorXYZ( level_, iterationVertex0, iterationVertex1, iterationVertex2 );
      for ( const auto & faceIdx : indexing::FaceIterator( levelinfo::num_microedges_per_edge( level_ ) - 1 ))
      {
        auto cellIdx2 = *cellIterator2;

        faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X, localCellID )] =
        cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationX )];
        faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y, localCellID )] =
        cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationY )];
        faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY, localCellID )] =
        cellData[edgedof::macrocell::index( level_, cellIdx2.x(), cellIdx2.y(), cellIdx2.z(), dstEdgeOrientationXY )];

        cellIterator2++;

        auto cellIdx3 = *cellIterator3;

        faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XYZ, localCellID )] =
        cellData[edgedof::macrocell::index(
        level_, cellIdx3.x(), cellIdx3.y(), cellIdx3.z(), edgedof::EdgeDoFOrientation::XYZ )];

        cellIterator3++;
      }
    }

    WALBERLA_ASSERT( cellIterator == cellIterator.end());
  }
  this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Cell to Face" );
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packVertexForCell( const Vertex*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   // nothing to do here
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackCellFromVertex( Cell*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   // nothing to do here
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalVertexToCell( const Vertex*, Cell* ) const
{
   // nothing to do here
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packEdgeForCell( const Edge*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );
   const ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );

   WALBERLA_UNUSED( receiver );
   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver ) );

   for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
   {
      buffer << edgeData[edgedof::macroedge::index( level_, it.x() )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = receiver->getLocalEdgeID( sender );
   const uint_t iterationVertex0 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( {iterationVertex0, iterationVertex1} ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); i++ )
   {
      ValueType tmp;
      buffer >> tmp;
      cellData[edgedof::macrocell::index( level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )] =
          tmp;
      cellIterator++;
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver->getID() ) );

   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = receiver->getLocalEdgeID( sender->getID() );
   const uint_t iterationVertex0 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( {iterationVertex0, iterationVertex1} ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   const ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );

   for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
   {
      cellData[edgedof::macrocell::index( level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )] =
          edgeData[edgedof::macroedge::index( level_, it.x() )];
      cellIterator++;
   }
}

template class EdgeDoFPackInfo< double >;
template class EdgeDoFPackInfo< float >;
template class EdgeDoFPackInfo< int >;
template class EdgeDoFPackInfo< long >;
template class EdgeDoFPackInfo< long long >;

} // namespace hyteg
