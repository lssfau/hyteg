/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/generatedKernels/communicate_directly_vertexdof_cell_to_face.hpp"
#include "hyteg/p1functionspace/generatedKernels/communicate_directly_vertexdof_face_to_cell.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

namespace hyteg {
namespace volumedofspace {

template < typename ValueType >
class VolumeDoFPackInfo : public communication::PackInfo
{
 public:
   VolumeDoFPackInfo( std::shared_ptr< PrimitiveStorage >                                      storage,
                      uint_t                                                                   level,
                      const std::map< PrimitiveID, uint_t >&                                   numScalarsPerPrimitive,
                      PrimitiveDataID< FunctionMemory< ValueType >, Face >                     faceInnerDataID,
                      PrimitiveDataID< FunctionMemory< ValueType >, Cell >                     cellInnerDataID,
                      std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs,
                      std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs,
                      indexing::VolumeDoFMemoryLayout                                          memoryLayout )
   : storage_( storage )
   , level_( level )
   , numScalarsPerPrimitive_( numScalarsPerPrimitive )
   , faceInnerDataID_( faceInnerDataID )
   , cellInnerDataID_( cellInnerDataID )
   , faceGhostLayerDataIDs_( faceGhostLayerDataIDs )
   , cellGhostLayerDataIDs_( cellGhostLayerDataIDs )
   , memoryLayout_( memoryLayout )
   {}

   void packVertexForEdge( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const override;

   void packEdgeForVertex( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackVertexFromEdge( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const override;

   void packEdgeForFace( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackFaceFromEdge( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const override;

   void packFaceForEdge( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackEdgeFromFace( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const override;

   void packFaceForCell( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromFace( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const override;

   void packCellForFace( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackFaceFromCell( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalCellToFace( const Cell* sender, Face* receiver ) const override;

   void packVertexForCell( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromVertex( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const override;

   void packEdgeForCell( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromEdge( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const override;

   /// @name Face to Face
   ///@{
   void packFaceForFace( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const;

   void unpackFaceFromFace( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const;

   void communicateLocalFaceToFace( const Face* sender, Face* receiver ) const;
   ///@}

   /// @name Cell to Cell
   ///@{
   void packCellForCell( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const;

   void unpackCellFromCell( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const;

   void communicateLocalCellToCell( const Cell* sender, Cell* receiver ) const;
   ///@}

 private:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              level_;
   std::map< PrimitiveID, uint_t >     numScalarsPerPrimitive_;

   /// Each data ID stores the inner scalar fields for all levels.
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceInnerDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellInnerDataID_;

   /// One data ID per ghost-layer. Should be up to 3 in 2D and 4 in 3D.
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs_;
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs_;

   indexing::VolumeDoFMemoryLayout memoryLayout_;
};

/// @name Face to Face
///@{
template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packFaceForFace( const Face*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "Macro-face to macro-face packing not implemented!" );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackFaceFromFace( Face*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Macro-face to macro-face packing not implemented!" );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalFaceToFace( const Face* sender, Face* receiver ) const
{
   this->storage_->getTimingTree()->start( "DG - Face to Face" );

   WALBERLA_CHECK_GREATER(
       numScalarsPerPrimitive_.count( sender->getID() ), 0, "Don't know how many scalars there are to send per volume." );

   // Which local edges do we iterate on?

   uint_t senderLocalEdgeID   = std::numeric_limits< uint_t >::max();
   uint_t receiverLocalEdgeID = std::numeric_limits< uint_t >::max();

   for ( const auto& [edgeID, npid] : sender->getIndirectNeighborFaceIDsOverEdges() )
   {
      if ( receiver->getID() == npid )
      {
         senderLocalEdgeID = edgeID;
         break;
      }
   }

   WALBERLA_ASSERT_LESS_EQUAL( senderLocalEdgeID, 2, "Couldn't find receiver face in neighborhood." );

   for ( const auto& [edgeID, npid] : receiver->getIndirectNeighborFaceIDsOverEdges() )
   {
      if ( sender->getID() == npid )
      {
         receiverLocalEdgeID = edgeID;
         break;
      }
   }

   WALBERLA_ASSERT_LESS_EQUAL( receiverLocalEdgeID, 2, "Couldn't find sender face in neighborhood." );

   // Do we need to invert the iteration direction?
   // We simply check the orientation of the interface macro-edge from both sides.
   // If the orientation does not change or if it changes twice, we do not need to switch up the iteration.

   const auto senderEdgeOrienatation   = sender->getEdgeOrientation()[senderLocalEdgeID];
   const auto receiverEdgeOrienatation = receiver->getEdgeOrientation()[receiverLocalEdgeID];
   const auto finalSenderOrientation   = senderEdgeOrienatation * receiverEdgeOrienatation;

   const auto numMicroVolumes = levelinfo::num_microedges_per_edge( level_ );

   const ValueType* faceData = sender->getData( faceInnerDataID_ )->getPointer( level_ );
   ValueType*       glData   = receiver->getData( faceGhostLayerDataIDs_.at( receiverLocalEdgeID ) )->getPointer( level_ );

   const auto faceBoundaryDirection = hyteg::indexing::getFaceBoundaryDirection( senderLocalEdgeID, finalSenderOrientation );
   const auto ndofs                 = numScalarsPerPrimitive_.at( sender->getID() );

   uint_t glMicroVolumeIdx = 0;

   for ( const auto& it : hyteg::indexing::FaceBoundaryIterator( numMicroVolumes, faceBoundaryDirection, 0 ) )
   {
      for ( uint_t dof = 0; dof < ndofs; dof++ )
      {
         const auto senderIdx =
             volumedofspace::indexing::index( it.x(), it.y(), facedof::FaceType::GRAY, dof, ndofs, level_, memoryLayout_ );
         const auto senderVal = faceData[senderIdx];

         const auto receiverIdx = indexGhostLayerDirectly( glMicroVolumeIdx, dof, ndofs, level_, memoryLayout_ );

         glData[receiverIdx] = senderVal;
      }
      glMicroVolumeIdx++;
   }
   this->storage_->getTimingTree()->stop( "DG - Face to Face" );
}

///@}

/// @name Cell to Cell
///@{
template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packCellForCell( const Cell*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "Macro-cell to macro-cell packing not implemented!" );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackCellFromCell( Cell*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Macro-cell to macro-cell packing not implemented!" );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalCellToCell( const Cell* sender, Cell* receiver ) const
{
   this->storage_->getTimingTree()->start( "DG - Cell to Cell" );

   WALBERLA_CHECK_GREATER(
       numScalarsPerPrimitive_.count( sender->getID() ), 0, "Don't know how many scalars there are to send per volume." );

   // Which local faces do we iterate on?

   uint_t senderLocalFaceID   = std::numeric_limits< uint_t >::max();
   uint_t receiverLocalFaceID = std::numeric_limits< uint_t >::max();

   for ( const auto& [faceID, npid] : sender->getIndirectNeighborCellIDsOverFaces() )
   {
      if ( receiver->getID() == npid )
      {
         senderLocalFaceID = faceID;
         break;
      }
   }

   WALBERLA_ASSERT_LESS_EQUAL( senderLocalFaceID, 3, "Couldn't find receiver cell in neighborhood." );

   for ( const auto& [faceID, npid] : receiver->getIndirectNeighborCellIDsOverFaces() )
   {
      if ( sender->getID() == npid )
      {
         receiverLocalFaceID = faceID;
         break;
      }
   }

   WALBERLA_ASSERT_LESS_EQUAL( receiverLocalFaceID, 3, "Couldn't find sender cell in neighborhood." );

   // We iterate in standard direction on the sender side.
   // This means iteration from the smallest to second-smallest index in the inner loop, to the largest in the outer loop.

   // On the receiver side we need to find the corresponding local vertex IDs.
   // The iteration is performed with the cell boundary iterator. On the receiver side, we need to write to the ghost-layer,
   // which we can achieve via extended access by setting one of the logical cell indices to -1.
   // That index depends on the local face ID on the receiver side.

   auto senderLocalVertexIDsSet = hyteg::indexing::cellLocalFaceIDsToSpanningVertexIDs.at( senderLocalFaceID );

   std::vector< uint_t >      senderLocalVertexIDs;
   std::vector< PrimitiveID > vertexPIDs;
   for ( auto slvid : senderLocalVertexIDsSet )
   {
      senderLocalVertexIDs.push_back( slvid );
      vertexPIDs.push_back( sender->neighborVertices().at( slvid ) );
   }

   // Sorting the receiver local vertex IDs correctly by matching the vertex primitive IDs.
   std::vector< uint_t > receiverLocalVertexIDs;
   for ( auto vpid : vertexPIDs )
   {
      receiverLocalVertexIDs.push_back( receiver->getLocalVertexID( vpid ) );
   }

   // We need one iterator for each cell type at the boundary.
   // This first one always iterates of the WHITE_UP type, the second cell type depends on the local macro-face ID.
   auto cellIteratorSenderCellType0 = hyteg::indexing::CellBoundaryIterator(
       levelinfo::num_microedges_per_edge( level_ ), senderLocalVertexIDs[0], senderLocalVertexIDs[1], senderLocalVertexIDs[2] );
   auto cellIteratorSenderCellType1 = hyteg::indexing::CellBoundaryIterator( levelinfo::num_microedges_per_edge( level_ ) - 1,
                                                                             senderLocalVertexIDs[0],
                                                                             senderLocalVertexIDs[1],
                                                                             senderLocalVertexIDs[2] );

   auto cellIteratorReceiverCellType0 = hyteg::indexing::CellBoundaryIterator( levelinfo::num_microedges_per_edge( level_ ),
                                                                               receiverLocalVertexIDs[0],
                                                                               receiverLocalVertexIDs[1],
                                                                               receiverLocalVertexIDs[2] );
   auto cellIteratorReceiverCellType1 = hyteg::indexing::CellBoundaryIterator( levelinfo::num_microedges_per_edge( level_ ) - 1,
                                                                               receiverLocalVertexIDs[0],
                                                                               receiverLocalVertexIDs[1],
                                                                               receiverLocalVertexIDs[2] );

   // What remains is to find out the cell type for the second iterator.
   celldof::CellType otherCellType = celldof::CellType::WHITE_DOWN; // makes it crash by default, WHITE_DOWN is never at boundary
   switch ( senderLocalFaceID )
   {
   case 0:
      otherCellType = celldof::CellType::BLUE_UP;
      break;
   case 1:
      otherCellType = celldof::CellType::GREEN_UP;
      break;
   case 2:
      otherCellType = celldof::CellType::BLUE_DOWN;
      break;
   case 3:
      otherCellType = celldof::CellType::GREEN_DOWN;
      break;
   default:
      WALBERLA_ABORT( "Invalid local face ID." );
   }

   const ValueType* cellData = sender->getData( cellInnerDataID_ )->getPointer( level_ );
   ValueType*       glData   = receiver->getData( cellGhostLayerDataIDs_.at( receiverLocalFaceID ) )->getPointer( level_ );

   const auto ndofs = numScalarsPerPrimitive_.at( sender->getID() );

   // Iterating over WHITE_UP cells ...
   while ( cellIteratorSenderCellType0 != cellIteratorSenderCellType0.end() )
   {
      for ( uint_t dof = 0; dof < ndofs; dof++ )
      {
         const auto senderIdx = volumedofspace::indexing::index( cellIteratorSenderCellType0->x(),
                                                                 cellIteratorSenderCellType0->y(),
                                                                 cellIteratorSenderCellType0->z(),
                                                                 celldof::CellType::WHITE_UP,
                                                                 dof,
                                                                 ndofs,
                                                                 level_,
                                                                 memoryLayout_ );
         const auto senderVal = cellData[senderIdx];

         const auto receiverIdx = indexing::indexNeighborInGhostLayer( receiverLocalFaceID,
                                                                       cellIteratorReceiverCellType0->x(),
                                                                       cellIteratorReceiverCellType0->y(),
                                                                       cellIteratorReceiverCellType0->z(),
                                                                       celldof::CellType::WHITE_UP,
                                                                       dof,
                                                                       ndofs,
                                                                       level_,
                                                                       memoryLayout_ );

         glData[receiverIdx] = senderVal;
      }
      cellIteratorSenderCellType0++;
      cellIteratorReceiverCellType0++;
   }

   // ... and the other type.
   while ( cellIteratorSenderCellType1 != cellIteratorSenderCellType1.end() )
   {
      for ( uint_t dof = 0; dof < ndofs; dof++ )
      {
         const auto senderIdx = volumedofspace::indexing::index( cellIteratorSenderCellType1->x(),
                                                                 cellIteratorSenderCellType1->y(),
                                                                 cellIteratorSenderCellType1->z(),
                                                                 otherCellType,
                                                                 dof,
                                                                 ndofs,
                                                                 level_,
                                                                 memoryLayout_ );
         const auto senderVal = cellData[senderIdx];

         // Cheat!!!
         // The cell type is in general not identical on the receiver side! But we don't care as the ghost-layer indexing
         // either accepts WHITE_UP for the first type of cells, or any other cell type for the second type of cells.
         const auto receiverIdx = indexing::indexNeighborInGhostLayer( receiverLocalFaceID,
                                                                       cellIteratorReceiverCellType1->x(),
                                                                       cellIteratorReceiverCellType1->y(),
                                                                       cellIteratorReceiverCellType1->z(),
                                                                       otherCellType,
                                                                       dof,
                                                                       ndofs,
                                                                       level_,
                                                                       memoryLayout_ );

         glData[receiverIdx] = senderVal;
      }
      cellIteratorSenderCellType1++;
      cellIteratorReceiverCellType1++;
   }

   this->storage_->getTimingTree()->stop( "DG - Cell to Cell" );
}
///@}

/// @name Vertex to Edge
///@{
template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                        const PrimitiveID&         receiver,
                                                        walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                           const PrimitiveID&         sender,
                                                           walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

///@}
/// @name Edge to Vertex
///@{

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                        const PrimitiveID&         receiver,
                                                        walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                           const PrimitiveID&         sender,
                                                           walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

///@}
/// @name Edge to Face
///@{

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

///@}
/// @name Face to Edge
///@{

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template <>
inline void VolumeDoFPackInfo< real_t >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
inline void VolumeDoFPackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template <>
inline void VolumeDoFPackInfo< real_t >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
inline void VolumeDoFPackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packVertexForCell( const Vertex*              sender,
                                                        const PrimitiveID&         receiver,
                                                        walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackCellFromVertex( Cell*                      receiver,
                                                           const PrimitiveID&         sender,
                                                           walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::packEdgeForCell( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
   WALBERLA_UNUSED( buffer );
}

template < typename ValueType >
void VolumeDoFPackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   WALBERLA_ABORT(
       "This communication pattern is not meaningful for volume DoFs - should not have been called in the first place." );
   WALBERLA_UNUSED( sender );
   WALBERLA_UNUSED( receiver );
}

} // namespace volumedofspace
} //namespace hyteg
