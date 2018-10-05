#pragma once

#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"

namespace hhg {

using walberla::uint_t;

template < typename ValueType >
class EdgeDoFPackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   EdgeDoFPackInfo( uint_t                                                 level,
                    PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                    PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                    PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                    PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                    std::weak_ptr< PrimitiveStorage >                      storage )
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

   void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const override;

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDCell_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;
};

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   ValueType* edgeData       = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t     vertexIdOnEdge = sender->vertex_index( receiver );
   if( vertexIdOnEdge == 0 )
   {
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_HO_C )];
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_VE_SE )];
      if( sender->getNumNeighborFaces() == 2 )
      {
         buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_DI_N )];
      }
   } else if( vertexIdOnEdge == 1 )
   {
      uint_t nbrEdges = levelinfo::num_microedges_per_edge( level_ );
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_HO_C )];
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_DI_S )];
      if( sender->getNumNeighborFaces() == 2 )
      {
         buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_VE_NW )];
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
      vertexData[edgeIdOnVertex] =
          edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_HO_C )];
      vertexData[receiver->getNumNeighborEdges() + receiver->face_index( sender->neighborFaces()[0] )] =
          edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_VE_SE )];
      if( sender->getNumNeighborFaces() == 2 )
      {
         vertexData[receiver->getNumNeighborEdges() + receiver->face_index( sender->neighborFaces()[1] )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, 0, stencilDirection::EDGE_DI_N )];
      }
   } else if( vertexIdOnEdge == 1 )
   {
      uint_t nbrEdges = levelinfo::num_microedges_per_edge( level_ );
      vertexData[edgeIdOnVertex] =
          edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_HO_C )];
      vertexData[receiver->getNumNeighborEdges() + receiver->face_index( sender->neighborFaces()[0] )] =
          edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_DI_S )];
      if( sender->getNumNeighborFaces() == 2 )
      {
         vertexData[receiver->getNumNeighborEdges() + receiver->face_index( sender->neighborFaces()[1] )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, nbrEdges - 1, stencilDirection::EDGE_VE_NW )];
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
   ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
   for( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); ++i )
   {
      buffer << edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, i, stencilDirection::EDGE_HO_C )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   using edgedof::macroface::indexFromHorizontalEdge;
   using hhg::edgedof::macroface::BorderIterator;
   ValueType*                    faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = receiver->edge_index( sender );
   indexing::FaceBorderDirection faceDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
   for( const auto& it : BorderIterator( level_, faceDir, 0 ) )
   {
      if( edgeIndexOnFace == 0 )
      {
         buffer >>
             faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 1 )
      {
         buffer >>
             faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_DI_N )];
      } else if( edgeIndexOnFace == 2 )
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
   using edgedof::macroface::indexFromHorizontalEdge;
   using hhg::edgedof::macroface::BorderIterator;
   ValueType*                    faceData        = receiver->getData( dataIDFace_ )->getPointer( level_ );
   ValueType*                    edgeData        = sender->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = receiver->edge_index( sender->getID() );
   indexing::FaceBorderDirection faceDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
   uint_t indexOnEdge = 0;
   for( const auto& it : BorderIterator( level_, faceDir, 0 ) )
   {
      if( edgeIndexOnFace == 0 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_HO_C )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, indexOnEdge, stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 1 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_DI_N )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, indexOnEdge, stencilDirection::EDGE_HO_C )];
      } else if( edgeIndexOnFace == 2 )
      {
         faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), stencilDirection::EDGE_VE_NW )] =
             edgeData[edgedof::macroedge::indexFromHorizontalEdge( level_, indexOnEdge, stencilDirection::EDGE_HO_C )];
      } else
      {
         WALBERLA_ABORT( "Wrong edgeIndexOnFace" )
      }
      ++indexOnEdge;
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   using hhg::edgedof::macroface::BorderIterator;
   ValueType*                    faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = sender->edge_index( receiver );
   indexing::FaceBorderDirection faceBorderDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );
   stencilDirection faceDirOne;
   stencilDirection faceDirTwo;
   stencilDirection faceDirThree;
   if( edgeIndexOnFace == 0 )
   {
      faceDirOne = stencilDirection::EDGE_HO_C;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_DI_N;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_DI_N;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      }
   } else if( edgeIndexOnFace == 1 )
   {
      faceDirOne = stencilDirection::EDGE_DI_N;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_HO_C;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_HO_C;
      }
   } else if( edgeIndexOnFace == 2 )
   {
      faceDirOne = stencilDirection::EDGE_VE_NW;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
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
   for( const auto& it : BorderIterator( level_, faceBorderDir, 1 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirOne )];
   }
   for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirTwo )];
   }
   for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
   {
      buffer << faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirThree )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType*       edgeData      = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t           faceIdOnEdge  = receiver->face_index( sender );
   stencilDirection dirHorizontal = faceIdOnEdge == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;
   /// first edge is south edge by convention
   for( uint_t i = 1; i < levelinfo::num_microvertices_per_edge( level_ ) - 1; ++i )
   {
      buffer >> edgeData[edgedof::macroedge::indexFromVertex( level_, i, dirHorizontal )];
   }
   stencilDirection dirDiagonal = faceIdOnEdge == 0 ? stencilDirection::EDGE_DI_SW : stencilDirection::EDGE_VE_NW;
   for( uint_t i = 1; i < levelinfo::num_microvertices_per_edge( level_ ); ++i )
   {
      buffer >> edgeData[edgedof::macroedge::indexFromVertex( level_, i, dirDiagonal )];
   }
   stencilDirection dirVertical = faceIdOnEdge == 0 ? stencilDirection::EDGE_VE_S : stencilDirection::EDGE_DI_NW;
   for( uint_t i = 1; i < levelinfo::num_microvertices_per_edge( level_ ); ++i )
   {
      buffer >> edgeData[edgedof::macroedge::indexFromVertex( level_, i, dirVertical )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   using hhg::edgedof::macroface::BorderIterator;
   ValueType*                    faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = sender->edge_index( receiver->getID() );
   indexing::FaceBorderDirection faceBorderDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );
   stencilDirection faceDirOne;
   stencilDirection faceDirTwo;
   stencilDirection faceDirThree;
   if( edgeIndexOnFace == 0 )
   {
      faceDirOne = stencilDirection::EDGE_HO_C;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_DI_N;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_DI_N;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      }
   } else if( edgeIndexOnFace == 1 )
   {
      faceDirOne = stencilDirection::EDGE_DI_N;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
      {
         faceDirTwo   = stencilDirection::EDGE_HO_C;
         faceDirThree = stencilDirection::EDGE_VE_NW;
      } else
      {
         faceDirTwo   = stencilDirection::EDGE_VE_NW;
         faceDirThree = stencilDirection::EDGE_HO_C;
      }
   } else if( edgeIndexOnFace == 2 )
   {
      faceDirOne = stencilDirection::EDGE_VE_NW;
      if( sender->edge_orientation[edgeIndexOnFace] == 1 )
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

   ValueType*       edgeData      = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   uint_t           faceIdOnEdge  = receiver->face_index( sender->getID() );
   stencilDirection dirHorizontal = faceIdOnEdge == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;

   uint_t indexOnEdge = 1;
   for( const auto& it : BorderIterator( level_, faceBorderDir, 1 ) )
   {
      edgeData[edgedof::macroedge::indexFromVertex( level_, indexOnEdge, dirHorizontal )] =
          faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirOne )];
      ++indexOnEdge;
   }
   stencilDirection dirDiagonal = faceIdOnEdge == 0 ? stencilDirection::EDGE_DI_SW : stencilDirection::EDGE_VE_NW;
   indexOnEdge                  = 1;
   for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
   {
      edgeData[edgedof::macroedge::indexFromVertex( level_, indexOnEdge, dirDiagonal )] =
          faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirTwo )];
      ++indexOnEdge;
   }
   stencilDirection dirVertical = faceIdOnEdge == 0 ? stencilDirection::EDGE_VE_S : stencilDirection::EDGE_DI_NW;
   indexOnEdge                  = 1;
   for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
   {
      edgeData[edgedof::macroedge::indexFromVertex( level_, indexOnEdge, dirVertical )] =
          faceData[edgedof::macroface::indexFromHorizontalEdge( level_, it.col(), it.row(), faceDirThree )];
      ++indexOnEdge;
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   const ValueType* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
   ValueType*       cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = receiver->getLocalFaceID( sender->getID() );
   const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );
   auto         dstEdgeOrientationX =
       edgedof::convertEdgeDoFOrientation( edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationY =
       edgedof::convertEdgeDoFOrientation( edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientation(
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
}

} //namespace hhg
