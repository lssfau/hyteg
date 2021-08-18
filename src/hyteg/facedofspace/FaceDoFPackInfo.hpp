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

#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/facedofspace/FaceDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {

using facedof::macroface::CELL_BLUE;
using facedof::macroface::CELL_GRAY;
using walberla::uint_t;

template < typename ValueType >
class FaceDoFPackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   FaceDoFPackInfo( uint_t                                                 level,
                    PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                    PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                    PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                    std::weak_ptr< PrimitiveStorage >                      storage )
   : communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, storage )
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

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;
};

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   /// see DGMemory.hpp for a description of the Vertex Memory
   ValueType* vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );
   for ( const PrimitiveID& faceID : storage_.lock()->getEdge( receiver )->neighborFaces() )
   {
      buffer << vertexData[sender->face_index( faceID ) * 2];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   typedef stencilDirection sD;
   ValueType*               edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t                   pos      = std::numeric_limits< uint_t >::max();
   if ( receiver->vertex_index( sender ) == 0 )
   {
      pos = 0;
   }
   else if ( receiver->vertex_index( sender ) == 1 )
   {
      pos = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   else
   {
      WALBERLA_LOG_WARNING( "Vertex with ID: " << sender.getID() << " is not in Edge: " << receiver )
   }
   buffer >> edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_GRAY_SE )];
   if ( receiver->getNumNeighborFaces() == 2 )
   {
      buffer >> edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_GRAY_NE )];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   typedef stencilDirection sD;
   ValueType*               vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );
   ValueType*               edgeData   = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t                   pos        = std::numeric_limits< uint_t >::max();
   if ( receiver->vertex_index( sender->getID() ) == 0 )
   {
      pos = 0;
   }
   else if ( receiver->vertex_index( sender->getID() ) == 1 )
   {
      pos = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   else
   {
      WALBERLA_LOG_WARNING( "Vertex with ID: " << sender << " is not in Edge: " << receiver )
   }
   edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_GRAY_SE )] =
       vertexData[sender->face_index( receiver->neighborFaces()[0] ) * 2];
   if ( receiver->getNumNeighborFaces() == 2 )
   {
      edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_GRAY_NE )] =
          vertexData[sender->face_index( receiver->neighborFaces()[1] ) * 2];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   ///the blue face DoF which are owned by the face need to communicated to the vertex
   typedef stencilDirection sD;
   ValueType*               edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t                   pos      = std::numeric_limits< uint_t >::max();
   if ( sender->vertex_index( receiver ) == 0 )
   {
      pos = 1;
   }
   else if ( sender->vertex_index( receiver ) == 1 )
   {
      pos = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   else
   {
      WALBERLA_LOG_WARNING( "Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender )
   }
   buffer << edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_BLUE_SE )];
   if ( sender->getNumNeighborFaces() == 2 )
   {
      buffer << edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, sD::CELL_BLUE_NW )];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   for ( const PrimitiveID& faceID : storage_.lock()->getEdge( sender )->neighborFaces() )
   {
      buffer >> vertexData[receiver->face_index( faceID ) * 2 + 1];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     pos      = std::numeric_limits< uint_t >::max();
   if ( sender->vertex_index( receiver->getID() ) == 0 )
   {
      pos = 1;
   }
   else if ( sender->vertex_index( receiver->getID() ) == 1 )
   {
      pos = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   else
   {
      WALBERLA_LOG_WARNING( "Vertex with ID: " << receiver << " is not in Edge: " << sender )
   }
   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   vertexData[receiver->face_index( sender->neighborFaces()[0] ) * 2 + 1] =
       edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, stencilDirection::CELL_BLUE_SE )];
   ;
   if ( sender->getNumNeighborFaces() == 2 )
   {
      vertexData[receiver->face_index( sender->neighborFaces()[1] ) * 2 + 1] =
          edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, stencilDirection::CELL_BLUE_NW )];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   ValueType*       edgeData     = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t           vPerEdge     = levelinfo::num_microvertices_per_edge( level_ );
   uint_t           faceIdOnEdge = sender->face_index( receiver );
   stencilDirection dirCellGray;
   //the first face is the south face and the second the north face
   if ( faceIdOnEdge == 0 )
   {
      dirCellGray = stencilDirection::CELL_GRAY_SE;
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( faceIdOnEdge, 1 );
      dirCellGray = stencilDirection::CELL_GRAY_NE;
   }
   for ( uint_t i = 0; i < vPerEdge - 1; ++i )
   {
      buffer << edgeData[facedof::macroedge::indexFaceFromVertex( level_, i, dirCellGray )];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType* faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t     edgeIndexOnFace = receiver->edge_index( sender );
   for ( auto it =
             facedof::macroface::indexIterator( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_ );
         it != facedof::macroface::indexIterator();
         ++it )
   {
      buffer >> faceData[*it];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   ValueType*       edgeData        = sender->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType*       faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t           faceIdOnEdge    = sender->face_index( receiver->getID() );
   uint_t           edgeIndexOnFace = receiver->edge_index( sender->getID() );
   stencilDirection dirCellGray;
   //the first face is the south face and the second the north face
   if ( faceIdOnEdge == 0 )
   {
      dirCellGray = stencilDirection::CELL_GRAY_SE;
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( faceIdOnEdge, 1 );
      dirCellGray = stencilDirection::CELL_GRAY_NE;
   }
   uint_t pos = 0;
   for ( auto it =
             facedof::macroface::indexIterator( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_ );
         it != facedof::macroface::indexIterator();
         ++it )
   {
      faceData[*it] = edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, dirCellGray )];
      pos++;
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   ValueType* faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t     edgeIndexOnFace = sender->edge_index( receiver );
   for ( auto it =
             facedof::macroface::indexIterator( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_ );
         it != facedof::macroface::indexIterator();
         ++it )
   {
      buffer << faceData[*it];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType*       edgeData     = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t           vPerEdge     = levelinfo::num_microvertices_per_edge( level_ );
   uint_t           faceIdOnEdge = receiver->face_index( sender );
   stencilDirection dirCellBlue;
   //the first face is the south face and the second the north face
   if ( faceIdOnEdge == 0 )
   {
      dirCellBlue = stencilDirection::CELL_BLUE_SE;
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( faceIdOnEdge, 1 );
      dirCellBlue = stencilDirection::CELL_BLUE_NW;
   }
   //unpack Blue Cell
   for ( uint_t i = 1; i < vPerEdge - 1; ++i )
   {
      buffer >> edgeData[facedof::macroedge::indexFaceFromVertex( level_, i, dirCellBlue )];
   }
}

template < typename ValueType >
void FaceDoFPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   ValueType*       edgeData        = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType*       faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   const uint_t     faceIdOnEdge    = receiver->face_index( sender->getID() );
   const uint_t     edgeIndexOnFace = sender->edge_index( receiver->getID() );
   uint_t           pos             = 1;
   stencilDirection dirCellBlue;
   //the first face is the south face and the second the north face
   if ( faceIdOnEdge == 0 )
   {
      dirCellBlue = stencilDirection::CELL_BLUE_SE;
   }
   else
   {
      WALBERLA_ASSERT_EQUAL( faceIdOnEdge, 1 );
      dirCellBlue = stencilDirection::CELL_BLUE_NW;
   }
   for ( auto it =
             facedof::macroface::indexIterator( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_ );
         it != facedof::macroface::indexIterator();
         ++it )
   {
      edgeData[facedof::macroedge::indexFaceFromVertex( level_, pos, dirCellBlue )] = faceData[*it];
      pos++;
   }
}

} //namespace hyteg
