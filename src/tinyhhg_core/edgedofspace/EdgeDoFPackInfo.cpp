#include <tinyhhg_core/Algorithms.hpp>
#include "EdgeDoFPackInfo.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"

namespace hhg {

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
   WALBERLA_UNUSED( receiver );
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

//   for ( uint_t neighborCellID = 0; neighborCellID < sender->getNumNeighborCells(); neighborCellID++  )
//   {
//      const Cell & neighborCell = *( storage_.lock()->getCell( sender->neighborCells().at( neighborCellID ) ) );

   for( const auto edgeOrientationOnReferenceEdge: { hhg::edgedof::EdgeDoFOrientation::Z,hhg::edgedof::EdgeDoFOrientation::YZ})
   {
      indexOnEdge = 0;
      for( const auto &neighborCellID : sender->neighborCells() )
      {
         const Cell& neighborCell    = *( storage_.lock()->getCell( neighborCellID ) );
         const auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( receiver->getID() );
         const auto edgeLocalCellID = receiver->cell_index( neighborCellID);
         const auto faceLocalCellID = sender->cell_index( neighborCellID);

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
              neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

         const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationCellToFace(
             edgeOrientationOnReferenceEdge, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

      for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
      {
         edgeData[edgedof::macroedge::indexOnNeighborCell( level_, indexOnEdge, edgeLocalCellID, receiver->getNumNeighborFaces(),edgeOrientationOnReferenceEdge  )] =
            faceData[edgedof::macroface::index( level_, it.col(), it.row(), edgeOrientationOnFace, faceLocalCellID )];
         ++indexOnEdge;
      }
      }
   }
   for( const auto edgeOrientationOnReferenceEdge: { hhg::edgedof::EdgeDoFOrientation::Z,hhg::edgedof::EdgeDoFOrientation::YZ})
   {
      indexOnEdge = 0;
      for( const auto &neighborCellID : sender->neighborCells() )
      {
         const Cell& neighborCell    = *( storage_.lock()->getCell( neighborCellID ) );
         const auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( receiver->getID() );
         const auto edgeLocalCellID = receiver->cell_index( neighborCellID);
         const auto faceLocalCellID = sender->cell_index( neighborCellID);

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
            {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
             neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

         const auto edgeOrientationOnFace = edgedof::convertEdgeDoFOrientationCellToFace(
            edgeOrientationOnReferenceEdge, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

         for( const auto& it : BorderIterator( level_, faceBorderDir, 0 ) )
         {
            edgeData[edgedof::macroedge::indexOnNeighborCell( level_, indexOnEdge, edgeLocalCellID, receiver->getNumNeighborFaces(),edgeOrientationOnReferenceEdge  )] =
               faceData[edgedof::macroface::index( level_, it.col(), it.row(), edgeOrientationOnFace, faceLocalCellID )];
            ++indexOnEdge;
         }
      }
   }

}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_UNUSED( receiver );
   const ValueType* faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
   for( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
      buffer << faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = receiver->getLocalFaceID( sender );
   const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationY =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2);

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

   auto dstEdgeOrientationX =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationY =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2);

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

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   const ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationY =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2);

   auto dstEdgeOrientationZ =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Z, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationYZ = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::YZ, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXZ = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XZ, iterationVertex0, iterationVertex1, iterationVertex2);

   for( const auto& cellIt :
      edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationZ )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationYZ )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationXZ )];
   }

   auto cellItXYZ = edgedof::macrocell::BoundaryIteratorXYZ( level_, iterationVertex0, iterationVertex1, iterationVertex2 );
   for( const auto& cellIt :
      edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 ) )
   {
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationX )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationY )];
      buffer << cellData[edgedof::macrocell::index( level_, cellIt.x(), cellIt.y(), cellIt.z(), dstEdgeOrientationXY )];

      buffer << cellData[edgedof::macrocell::index(
         level_, cellItXYZ->x(), cellItXYZ->y(), cellItXYZ->z(), edgedof::EdgeDoFOrientation::XYZ )];
      cellItXYZ++;
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
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

   for( const auto& faceIdx : indexing::FaceIterator( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
   {
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::X, localCellID )];
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::Y, localCellID )];
      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::XY, localCellID )];

      buffer >> faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), EO::XYZ, localCellID )];
   }
}

template < typename ValueType >
void EdgeDoFPackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   const ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
   ValueType*       faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver->getID() );
   const uint_t localCellID      = receiver->cell_index( sender->getID() );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto dstEdgeOrientationX =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationY =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2);

   auto dstEdgeOrientationZ =
      edgedof::convertEdgeDoFOrientationFaceToCell(edgedof::EdgeDoFOrientation::Z, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationYZ = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::YZ, iterationVertex0, iterationVertex1, iterationVertex2);
   auto dstEdgeOrientationXZ = edgedof::convertEdgeDoFOrientationFaceToCell(
      edgedof::EdgeDoFOrientation::XZ, iterationVertex0, iterationVertex1, iterationVertex2);

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

   auto cellIterator2 = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );
   auto cellIterator3 = edgedof::macrocell::BoundaryIteratorXYZ( level_, iterationVertex0, iterationVertex1, iterationVertex2 );
   for( const auto& faceIdx : indexing::FaceIterator( levelinfo::num_microedges_per_edge( level_ ) - 1 ) )
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

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
}

template class EdgeDoFPackInfo< double >;
template class EdgeDoFPackInfo< float >;
template class EdgeDoFPackInfo< int >;
template class EdgeDoFPackInfo< uint_t >;

}/// namespace hhg

