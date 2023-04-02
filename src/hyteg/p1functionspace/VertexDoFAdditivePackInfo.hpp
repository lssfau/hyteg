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

#include "hyteg/communication/DoFSpacePackInfo.hpp"
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
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

///@}
/// @name Edge to Vertex
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                                const PrimitiveID&         receiver,
                                                                walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                                   const PrimitiveID&         sender,
                                                                   walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

///@}
/// @name Edge to Face
///@{

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::packEdgeForFace( const Edge* sender,
                                                              const PrimitiveID& /*receiver*/,
                                                              walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                                 const PrimitiveID&         sender,
                                                                 walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template < typename ValueType >
void VertexDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
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
      buffer << faceData[vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C )];
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
          faceData[vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C )];
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
