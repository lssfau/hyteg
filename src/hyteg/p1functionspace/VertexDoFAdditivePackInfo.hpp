/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Nils Kohl, Andreas Burkhart.
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
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {

template < typename ValueType >
class VertexDoFAdditivePackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   VertexDoFAdditivePackInfo( uint_t                                                 level,
                              PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                              PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                              PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                              PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                              const std::weak_ptr< PrimitiveStorage >&               storage )
   : communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
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

   void packFaceForVertex( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackVertexFromFace( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const override;

   void packFaceForCell( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromFace( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const override;

   void packCellForFace( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackFaceFromCell( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalCellToFace( const Cell* sender, Face* receiver ) const override;

   void packCellForEdge( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackEdgeFromCell( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const override;

   void packCellForVertex( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackVertexFromCell( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalCellToVertex( const Cell* sender, Vertex* receiver ) const override;

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDCell_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;
};

/// @name Vertex to Edge
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                                const PrimitiveID&         receiver,
                                                                walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Vertex -> Edge only meaningful in 2D." );

   ValueType* vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );

   buffer << vertexData[0];
   for ( size_t i = 0; i < sender->getNumNeighborEdges(); i++ )
   {
      buffer << vertexData[i];
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Vertex -> Edge only meaningful in 2D." );

   Vertex&    vertex   = *( storage_.lock()->getVertex( sender ) );
   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType  tmp;

   // unpack buffer
   std::vector< ValueType > senderbuffer( vertex.getNumNeighborEdges() + 1 );
   for ( size_t i = 0; i < vertex.getNumNeighborEdges() + 1; i++ )
   {
      buffer >> senderbuffer[i];
   }

   // handle vertex on edge
   const uint_t vertexIdOnEdge = receiver->vertex_index( sender );
   uint_t       vertexIndex;
   uint_t       vertexIndexInner;
   if ( vertexIdOnEdge == 0 )
   {
      vertexIndex      = 0;
      vertexIndexInner = 1;
   }
   else
   {
      vertexIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      vertexIndexInner = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndex, stencilDirection::VERTEX_C )] += senderbuffer[0];

   // edge itself
   uint_t edgeIndexOnVertex = vertex.edge_index( receiver->getID() );
   edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndexInner, stencilDirection::VERTEX_C )] +=
       senderbuffer[edgeIndexOnVertex + 1];

   // iterate over edge neighbour faces
   for ( uint_t neighborFace = 0; neighborFace < receiver->getNumNeighborFaces(); neighborFace++ )
   {
      Face& face = *( storage_.lock()->getFace( receiver->neighborFaces()[neighborFace] ) );
      for ( auto& e : face.neighborEdges() )
      {
         Edge* neighbourEdgePtr = storage_.lock()->getEdge( e );
         // in the case of multiple ranks this just ignores edges we do not need (ones with no connection to the vertex)
         if ( neighbourEdgePtr != nullptr )
         {
            Edge& neighbourEdge = *( neighbourEdgePtr );

            if ( ( neighbourEdge.getID() != receiver->getID() ) &&
                 ( ( neighbourEdge.getVertexID0() == sender ) || ( neighbourEdge.getVertexID1() == sender ) ) )
            {
               uint_t edgeIndexOnVertexLocal = vertex.edge_index( neighbourEdge.getID() );
               if ( vertexIdOnEdge == 0 )
               {
                  edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_W )] +=
                      senderbuffer[edgeIndexOnVertexLocal + 1];
               }
               else if ( vertexIdOnEdge == 1 )
               {
                  edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_E )] +=
                      senderbuffer[edgeIndexOnVertexLocal + 1];
               }
            }
         }
      }
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Vertex -> Edge only meaningful in 2D." );

   // #################################
   // ###### handle sending side ######
   // #################################
   ValueType* vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );

   // ###################################
   // ###### handle receiving side ######
   // ###################################
   Vertex&    vertex   = *( storage_.lock()->getVertex( sender->getID() ) );
   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   // handle vertex on edge
   const uint_t vertexIdOnEdge = receiver->vertex_index( sender->getID() );
   uint_t       vertexIndex;
   uint_t       vertexIndexInner;
   if ( vertexIdOnEdge == 0 )
   {
      vertexIndex      = 0;
      vertexIndexInner = 1;
   }
   else
   {
      vertexIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      vertexIndexInner = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndex, stencilDirection::VERTEX_C )] += vertexData[0];

   // edge itself
   uint_t edgeIndexOnVertex = vertex.edge_index( receiver->getID() );
   edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndexInner, stencilDirection::VERTEX_C )] +=
       vertexData[edgeIndexOnVertex + 1];

   // iterate over edge neighbour faces
   for ( uint_t neighborFace = 0; neighborFace < receiver->getNumNeighborFaces(); neighborFace++ )
   {
      Face& face = *( storage_.lock()->getFace( receiver->neighborFaces()[neighborFace] ) );
      for ( auto& e : face.neighborEdges() )
      {
         Edge* neighbourEdgePtr = storage_.lock()->getEdge( e );
         // in the case of multiple ranks this just ignores edges we do not need (ones with no connection to the vertex)
         if ( neighbourEdgePtr != nullptr )
         {
            Edge& neighbourEdge = *( neighbourEdgePtr );

            if ( ( neighbourEdge.getID() != receiver->getID() ) &&
                 ( ( neighbourEdge.getVertexID0() == sender->getID() ) || ( neighbourEdge.getVertexID1() == sender->getID() ) ) )
            {
               uint_t edgeIndexOnVertexLocal = vertex.edge_index( neighbourEdge.getID() );
               if ( vertexIdOnEdge == 0 )
               {
                  edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_W )] +=
                      vertexData[edgeIndexOnVertexLocal + 1];
               }
               else if ( vertexIdOnEdge == 1 )
               {
                  edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_E )] +=
                      vertexData[edgeIndexOnVertexLocal + 1];
               }
            }
         }
      }
   }
}

///@}
/// @name Edge to Vertex
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                                const PrimitiveID&         receiver,
                                                                walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Vertex only meaningful in 2D." );

   Vertex&    vertex   = *( storage_.lock()->getVertex( receiver ) );
   ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType  tmp;

   // create vector of correct size
   std::vector< ValueType > senderbuffer( vertex.getNumNeighborEdges() + 1 );

   // handle vertex on edge
   const uint_t vertexIdOnEdge = sender->vertex_index( receiver );
   uint_t       vertexIndex;
   uint_t       vertexIndexInner;
   if ( vertexIdOnEdge == 0 )
   {
      vertexIndex      = 0;
      vertexIndexInner = 1;
   }
   else
   {
      vertexIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      vertexIndexInner = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   senderbuffer[0] = edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndex, stencilDirection::VERTEX_C )];

   // edge itself
   uint_t edgeIndexOnVertex = vertex.edge_index( sender->getID() );
   senderbuffer[edgeIndexOnVertex + 1] =
       edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndexInner, stencilDirection::VERTEX_C )];

   // iterate over edge neighbour faces
   for ( uint_t neighborFace = 0; neighborFace < sender->getNumNeighborFaces(); neighborFace++ )
   {
      Face& face = *( storage_.lock()->getFace( sender->neighborFaces()[neighborFace] ) );
      for ( auto& e : face.neighborEdges() )
      {
         Edge* neighbourEdgePtr = storage_.lock()->getEdge( e );
         // in the case of multiple ranks this just ignores edges we do not need (ones with no connection to the vertex)
         if ( neighbourEdgePtr != nullptr )
         {
            Edge& neighbourEdge = *( neighbourEdgePtr );

            if ( ( neighbourEdge.getID() != sender->getID() ) &&
                 ( ( neighbourEdge.getVertexID0() == receiver ) || ( neighbourEdge.getVertexID1() == receiver ) ) )
            {
               uint_t edgeIndexOnVertexLocal = vertex.edge_index( neighbourEdge.getID() );
               if ( vertexIdOnEdge == 0 )
               {
                  senderbuffer[edgeIndexOnVertexLocal + 1] = edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_W )];
               }
               else if ( vertexIdOnEdge == 1 )
               {
                  senderbuffer[edgeIndexOnVertexLocal + 1] = edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_E )];
               }
            }
         }
      }
   }

   // pack buffer
   for ( size_t i = 0; i < vertex.getNumNeighborEdges() + 1; i++ )
   {
      buffer << senderbuffer[i];
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Vertex only meaningful in 2D." );

   ValueType* vertexData     = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   uint_t     edgeIdOnVertex = receiver->edge_index( sender );
   ValueType  tmp;

   for ( size_t i = 0; i < receiver->getNumNeighborEdges() + 1; i++ )
   {
      buffer >> tmp;
      vertexData[i] += tmp;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Vertex only meaningful in 2D." );

   ValueType* edgeData   = sender->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   ValueType  tmp;

   // create vector of correct size
   std::vector< ValueType > senderbuffer( receiver->getNumNeighborEdges() + 1 );

   // handle vertex on edge
   const uint_t vertexIdOnEdge = sender->vertex_index( receiver->getID() );
   uint_t       vertexIndex;
   uint_t       vertexIndexInner;
   if ( vertexIdOnEdge == 0 )
   {
      vertexIndex      = 0;
      vertexIndexInner = 1;
   }
   else
   {
      vertexIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      vertexIndexInner = levelinfo::num_microvertices_per_edge( level_ ) - 2;
   }
   vertexData[0] += edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndex, stencilDirection::VERTEX_C )];

   // edge itself
   uint_t edgeIndexOnVertex = receiver->edge_index( sender->getID() );
   vertexData[edgeIndexOnVertex + 1] +=
       edgeData[vertexdof::macroedge::indexFromVertex( level_, vertexIndexInner, stencilDirection::VERTEX_C )];

   // iterate over edge neighbour faces
   for ( uint_t neighborFace = 0; neighborFace < sender->getNumNeighborFaces(); neighborFace++ )
   {
      Face& face = *( storage_.lock()->getFace( sender->neighborFaces()[neighborFace] ) );
      for ( auto& e : face.neighborEdges() )
      {
         Edge* neighbourEdgePtr = storage_.lock()->getEdge( e );
         // in the case of multiple ranks this just ignores edges we do not need (ones with no connection to the vertex)
         if ( neighbourEdgePtr != nullptr )
         {
            Edge& neighbourEdge = *( neighbourEdgePtr );

            if ( ( neighbourEdge.getID() != sender->getID() ) && ( ( neighbourEdge.getVertexID0() == receiver->getID() ) ||
                                                                   ( neighbourEdge.getVertexID1() == receiver->getID() ) ) )
            {
               uint_t edgeIndexOnVertexLocal = receiver->edge_index( neighbourEdge.getID() );
               if ( vertexIdOnEdge == 0 )
               {
                  vertexData[edgeIndexOnVertexLocal + 1] += edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_W )];
               }
               else if ( vertexIdOnEdge == 1 )
               {
                  vertexData[edgeIndexOnVertexLocal + 1] += edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
                      level_, vertexIndexInner, neighborFace, stencilDirection::VERTEX_E )];
               }
            }
         }
      }
   }
}

///@}
/// @name Edge to Face
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                              const PrimitiveID&         receiver,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Face only meaningful in 2D." );

   ValueType* edgeData  = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     v_perEdge = levelinfo::num_microvertices_per_edge( level_ );

   // edge DoFs
   for ( uint_t i = 0; i < v_perEdge; i++ )
   {
      buffer << edgeData[vertexdof::macroedge::indexFromVertex( level_, i, stencilDirection::VERTEX_C )];
   }

   // inner face DoFs
   uint_t faceIndexOnEdge = sender->face_index( receiver );
   for ( uint_t i = 0; i < v_perEdge - 1; i++ )
   {
      buffer << edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
          level_, i, faceIndexOnEdge, stencilDirection::VERTEX_E )];
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Face only meaningful in 2D." );

   ValueType*                      faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                          edgeIndexOnFace = receiver->edge_index( sender );
   indexing::FaceBoundaryDirection faceBorderDirection =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );

   ValueType tmp;

   // edge DoFs
   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0 ) )
   {
      buffer >> tmp;
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] += tmp;
   }

   // inner face DoFs
   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 1 ) )
   {
      buffer >> tmp;
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] += tmp;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Edge -> Face only meaningful in 2D." );

   ValueType* edgeData        = sender->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType* faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t     edgeIndexOnFace = receiver->edge_index( sender->getID() );

   // edge DoFs
   uint_t                          idx = 0;
   indexing::FaceBoundaryDirection faceBorderDirection =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );
   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0 ) )
   {
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] += edgeData[idx];
      idx++;
   }

   // inner face DoFs
   uint_t faceIndexOnEdge = sender->face_index( receiver->getID() );
   uint_t v_perEdge       = levelinfo::num_microvertices_per_edge( level_ );

   idx = 0;
   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 1 ) )
   {
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] +=
          edgeData[vertexdof::macroedge::indexFromVertexOnNeighborFace(
              level_, idx, faceIndexOnEdge, stencilDirection::VERTEX_E )];
      idx++;
   }
}

///@}
/// @name Face to Edge
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                              const PrimitiveID&         receiver,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   ValueType*                      faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                          edgeIndexOnFace = sender->edge_index( receiver );
   indexing::FaceBoundaryDirection faceBorderDirection =
       indexing::getFaceBoundaryDirection( edgeIndexOnFace, sender->getEdgeOrientation()[edgeIndexOnFace] );

   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0, 1 ) )
   {
      buffer << faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )];
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   const auto storage  = storage_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( storage.get() );
   WALBERLA_CHECK( storage->faceExistsLocally( sender ) || storage->faceExistsInNeighborhood( sender ) );

   for ( const auto& it : vertexdof::macroedge::Iterator( level_, 1 ) )
   {
      ValueType tmp;
      buffer >> tmp;
      edgeData[vertexdof::macroedge::index( level_, it.x() )] += tmp;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   ValueType*                      edgeData     = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType*                      faceData     = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                          edgeIdOnFace = sender->edge_index( receiver->getID() );
   indexing::FaceBoundaryDirection faceBorderDirection =
       indexing::getFaceBoundaryDirection( edgeIdOnFace, sender->getEdgeOrientation()[edgeIdOnFace] );
   vertexdof::macroedge::Iterator edgeIterator( level_, 1 );
   for ( const auto& it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0, 1 ) )
   {
      edgeData[vertexdof::macroedge::index( level_, edgeIterator->x() )] +=
          faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )];
      edgeIterator++;
   }
}

///@}
/// @name Face to Vertex
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packFaceForVertex( const Face*                sender,
                                                                const PrimitiveID&         receiver,
                                                                walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   ValueType*   faceData      = sender->getData( dataIDFace_ )->getPointer( level_ );
   const uint_t maxIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
   const uint_t localVertexID = sender->vertex_index( receiver );
   switch ( localVertexID )
   {
   case 0:
      buffer << faceData[vertexdof::macroface::indexFromVertex( level_, 0, 0, stencilDirection::VERTEX_C )];
      break;
   case 1:
      buffer << faceData[vertexdof::macroface::indexFromVertex( level_, maxIndex, 0, stencilDirection::VERTEX_C )];
      break;
   case 2:
      buffer << faceData[vertexdof::macroface::indexFromVertex( level_, 0, maxIndex, stencilDirection::VERTEX_C )];
      break;
   default:
      WALBERLA_ABORT( "Invalid local vertex ID." );
      break;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackVertexFromFace( Vertex*                    receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   ValueType  tmp;
   buffer >> tmp;
   vertexData[0] += tmp;
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );

   ValueType*   faceData      = sender->getData( dataIDFace_ )->getPointer( level_ );
   ValueType*   vertexData    = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   const uint_t maxIndex      = levelinfo::num_microvertices_per_edge( level_ ) - 1;
   const uint_t localVertexID = sender->vertex_index( receiver->getID() );
   switch ( localVertexID )
   {
   case 0:
      vertexData[0] += faceData[vertexdof::macroface::indexFromVertex( level_, 0, 0, stencilDirection::VERTEX_C )];
      break;
   case 1:
      vertexData[0] += faceData[vertexdof::macroface::indexFromVertex( level_, maxIndex, 0, stencilDirection::VERTEX_C )];
      break;
   case 2:
      vertexData[0] += faceData[vertexdof::macroface::indexFromVertex( level_, 0, maxIndex, stencilDirection::VERTEX_C )];
      break;
   default:
      WALBERLA_ABORT( "Invalid local vertex ID." );
      break;
   }
}

///@}
/// @name Face to Cell
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                              const PrimitiveID&         receiver,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

///@}
/// @name Cell to Face
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                              const PrimitiveID&         receiver,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Face only meaningful in 3D." );

   const ValueType* cellData         = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t     localFaceID      = sender->getLocalFaceID( receiver );
   const uint_t     iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t     iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t     iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   for ( const auto& it :
         vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 ) )
   {
      buffer << cellData[vertexdof::macrocell::indexFromVertex( level_, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )];
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Face only meaningful in 3D." );

   ValueType* faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

   for ( const auto& it : vertexdof::macroface::Iterator( level_ ) )
   {
      ValueType tmp;
      buffer >> tmp;
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] += tmp;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   const ValueType* cellData         = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t     localFaceID      = sender->getLocalFaceID( receiver->getID() );
   const uint_t     iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t     iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t     iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   ValueType* faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

   auto cellIterator = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 );

   for ( const auto& it : vertexdof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;
      faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C )] +=
          cellData[vertexdof::macrocell::indexFromVertex(
              level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C )];
      cellIterator++;
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
}

///@}
/// @name Cell to Edge
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packCellForEdge( const Cell*                sender,
                                                              const PrimitiveID&         receiver,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   const ValueType*   cellData                  = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t       localEdgeID               = sender->getLocalEdgeID( receiver );
   const uint_t       iterationVertex0          = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t       iterationVertex1          = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   std::set< uint_t > possibleIterationVertices = { 0, 1, 2, 3 };
   possibleIterationVertices.erase( iterationVertex0 );
   possibleIterationVertices.erase( iterationVertex1 );
   const uint_t iterationVertex2 = *possibleIterationVertices.begin();

   const uint_t edgeSize = levelinfo::num_microvertices_per_edge( level_ );
   auto         it = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 );
   for ( uint_t i = 0; i < edgeSize; i++ )
   {
      buffer << cellData[vertexdof::macrocell::indexFromVertex( level_, it->x(), it->y(), it->z(), stencilDirection::VERTEX_C )];
      it++;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackEdgeFromCell( Edge*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

   for ( const auto& it : vertexdof::macroedge::Iterator( level_ ) )
   {
      ValueType tmp;
      buffer >> tmp;
      edgeData[vertexdof::macroedge::indexFromVertex( level_, it.x(), stencilDirection::VERTEX_C )] += tmp;
   }
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   ValueType*         edgeData                  = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   const ValueType*   cellData                  = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t       localEdgeID               = sender->getLocalEdgeID( receiver->getID() );
   const uint_t       iterationVertex0          = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t       iterationVertex1          = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   std::set< uint_t > possibleIterationVertices = { 0, 1, 2, 3 };
   possibleIterationVertices.erase( iterationVertex0 );
   possibleIterationVertices.erase( iterationVertex1 );
   const uint_t iterationVertex2 = *possibleIterationVertices.begin();

   const uint_t edgeSize = levelinfo::num_microvertices_per_edge( level_ );
   auto         it = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 );
   for ( uint_t i = 0; i < edgeSize; i++ )
   {
      edgeData[i] +=
          cellData[vertexdof::macrocell::indexFromVertex( level_, it->x(), it->y(), it->z(), stencilDirection::VERTEX_C )];
      it++;
   }
}

///@}
/// @name Cell to Vertex
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packCellForVertex( const Cell*                sender,
                                                                const PrimitiveID&         receiver,
                                                                walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );

   const ValueType* cellData      = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t     localVertexID = sender->getLocalVertexID( receiver );
   indexing::Index  microVertexIndexInMacroCell( 0, 0, 0 );
   switch ( localVertexID )
   {
   case 1:
      microVertexIndexInMacroCell.x() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   case 2:
      microVertexIndexInMacroCell.y() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   case 3:
      microVertexIndexInMacroCell.z() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   default:
      break;
   }
   buffer << cellData[vertexdof::macrocell::indexFromVertex( level_,
                                                             microVertexIndexInMacroCell.x(),
                                                             microVertexIndexInMacroCell.y(),
                                                             microVertexIndexInMacroCell.z(),
                                                             stencilDirection::VERTEX_C )];
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackVertexFromCell( Vertex*                    receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );

   ValueType* vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

   ValueType tmp;
   buffer >> tmp;
   vertexData[0] += tmp;
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalCellToVertex( const Cell* sender, Vertex* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   ValueType*       vertexData    = receiver->getData( dataIDVertex_ )->getPointer( level_ );
   const ValueType* cellData      = sender->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t     localVertexID = sender->getLocalVertexID( receiver->getID() );
   indexing::Index  microVertexIndexInMacroCell( 0, 0, 0 );
   switch ( localVertexID )
   {
   case 1:
      microVertexIndexInMacroCell.x() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   case 2:
      microVertexIndexInMacroCell.y() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   case 3:
      microVertexIndexInMacroCell.z() = levelinfo::num_microvertices_per_edge( level_ ) - 1;
      break;
   default:
      break;
   }
   vertexData[0] += cellData[vertexdof::macrocell::indexFromVertex( level_,
                                                                    microVertexIndexInMacroCell.x(),
                                                                    microVertexIndexInMacroCell.y(),
                                                                    microVertexIndexInMacroCell.z(),
                                                                    stencilDirection::VERTEX_C )];
}

///@}

} //namespace hyteg
