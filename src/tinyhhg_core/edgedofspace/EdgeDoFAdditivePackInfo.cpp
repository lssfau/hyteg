#include "EdgeDoFAdditivePackInfo.hpp"

#include "tinyhhg_core/Algorithms.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/indexing/DistanceCoordinateSystem.hpp"
#include "tinyhhg_core/indexing/LocalIDMappings.hpp"

namespace hhg {

/// @name Vertex to Edge
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packVertexForEdge(const Vertex *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackEdgeFromVertex(Edge *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *, Edge *) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

///@}
/// @name Edge to Vertex
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packEdgeForVertex(const Edge *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackVertexFromEdge(Vertex *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *, Vertex *) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

///@}
/// @name Edge to Face
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packEdgeForFace(const Edge *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackFaceFromEdge(Face *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *, Face *) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

///@}
/// @name Face to Edge
///@{

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                            const PrimitiveID&         receiver,
                                                            walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );
   using hhg::edgedof::macroface::BoundaryIterator;
   ValueType*                    faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = sender->edge_index( receiver );
   indexing::FaceBoundaryDirection faceBorderDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );

   edgedof::EdgeDoFOrientation orientation;
   switch ( edgeIndexOnFace )
   {
   case 0:
      orientation = edgedof::EdgeDoFOrientation::X;
      break;
   case 1:
      orientation = edgedof::EdgeDoFOrientation::XY;
      break;
   case 2:
      orientation = edgedof::EdgeDoFOrientation::Y;
      break;
   default:
      WALBERLA_ABORT( "Invalid orienation" );
   }

   for ( const auto& it : BoundaryIterator( level_, faceBorderDir, 0, 0 ) )
   {
      buffer << faceData[edgedof::macroface::index( level_, it.x(), it.y(), orientation )];
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                               const PrimitiveID&,
                                                               walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );

   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   for ( const auto& it : edgedof::macroedge::Iterator( level_, 0 ) )
   {
      ValueType tmp;
      buffer >> tmp;

      if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) != boundaryTypeToSkip_ )
      {
         edgeData[edgedof::macroedge::index( level_, it.x() )] += tmp;
      }
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );

   if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) == boundaryTypeToSkip_ )
      return;

   using hhg::edgedof::macroface::BoundaryIterator;
   ValueType*                    edgeData        = receiver->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType*                    faceData        = sender->getData( dataIDFace_ )->getPointer( level_ );
   uint_t                        edgeIndexOnFace = sender->edge_index( receiver->getID() );
   indexing::FaceBoundaryDirection faceBorderDir =
       indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );

   edgedof::EdgeDoFOrientation orientation;
   switch ( edgeIndexOnFace )
   {
   case 0:
      orientation = edgedof::EdgeDoFOrientation::X;
      break;
   case 1:
      orientation = edgedof::EdgeDoFOrientation::XY;
      break;
   case 2:
      orientation = edgedof::EdgeDoFOrientation::Y;
      break;
   default:
      WALBERLA_ABORT( "Invalid orienation" );
   }

   edgedof::macroedge::Iterator edgeIterator( level_, 0 );
   for ( const auto& it : BoundaryIterator( level_, faceBorderDir, 0, 0 ) )
   {
      edgeData[edgedof::macroedge::index( level_, edgeIterator->x() )] +=
          faceData[edgedof::macroface::index( level_, it.x(), it.y(), orientation )];
      edgeIterator++;
   }
}

///@}
/// @name Face to Vertex
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packFaceForVertex(const Face *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackVertexFromFace(Vertex *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalFaceToVertex(const Face *, Vertex *) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

///@}
/// @name Face to Cell
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packFaceForCell(const Face *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackCellFromFace(Cell *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalFaceToCell(const Face *, Cell *) const
{
   WALBERLA_ABORT( "Additive communication Face -> Cell not supported." );
}

///@}
/// @name Cell to Face
///@{

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                            const PrimitiveID&         receiver,
                                                            walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Face only meaningful in 3D." );

   ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto srcEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto srcEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto srcEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( const auto& cellIterator :
         edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      buffer << cellData[edgedof::macrocell::index(
          level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationX )];
      buffer << cellData[edgedof::macrocell::index(
          level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationY )];
      buffer << cellData[edgedof::macrocell::index(
          level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationXY )];
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                               const PrimitiveID&         sender,
                                                               walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Face only meaningful in 3D." );

   WALBERLA_UNUSED( sender );
   ValueType* faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );
   for ( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      ValueType tmpX, tmpY, tmpXY;
      buffer >> tmpX;
      buffer >> tmpY;
      buffer >> tmpXY;
      if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) != boundaryTypeToSkip_ )
      {
         faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )] += tmpX;
         faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )] += tmpY;
         faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )] += tmpXY;
      }
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 3D." );

   if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) == boundaryTypeToSkip_ )
   {
      return;
   }

   ValueType*       faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );
   const ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localFaceID      = sender->getLocalFaceID( receiver->getID() );
   const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
   const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
   const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

   auto srcEdgeOrientationX = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto srcEdgeOrientationY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::Y, iterationVertex0, iterationVertex1, iterationVertex2 );
   auto srcEdgeOrientationXY = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::XY, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )] +=
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationX )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )] +=
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationY )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )] +=
          cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationXY )];
      cellIterator++;
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
}

///@}
/// @name Cell to Edge
///@{

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packCellForEdge( const Cell*                sender,
                                                            const PrimitiveID&         receiver,
                                                            walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = sender->getLocalEdgeID( receiver );
   const uint_t iterationVertex0 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( {iterationVertex0, iterationVertex1} ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); i++ )
   {
      buffer << cellData[edgedof::macrocell::index(
          level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )];
      cellIterator++;
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackEdgeFromCell( Edge*                      receiver,
                                                               const PrimitiveID&         sender,
                                                               walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );
   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

   if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) == boundaryTypeToSkip_ )
   {
      for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
      {
         WALBERLA_UNUSED( it );
         ValueType tmp;
         buffer >> tmp;
      }
   }
   else
   {
      for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
      {
         ValueType tmp;
         buffer >> tmp;
         edgeData[edgedof::macroedge::index( level_, it.x() )] += tmp;
      }
   }
}

template < typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   if ( boundaryCondition_.getBoundaryType( receiver->getMeshBoundaryFlag() ) == boundaryTypeToSkip_ )
   {
      return;
   }

   ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = sender->getLocalEdgeID( receiver->getID() );
   const uint_t iterationVertex0 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( {iterationVertex0, iterationVertex1} ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
   {
      edgeData[edgedof::macroedge::index( level_, it.x() )] += cellData[edgedof::macrocell::index(
          level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )];
      cellIterator++;
   }
}

///@}
/// @name Cell to Vertex
///@{

template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::packCellForVertex(const Cell *, const PrimitiveID &, walberla::mpi::SendBuffer &) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}


template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::unpackVertexFromCell(Vertex *, const PrimitiveID &, walberla::mpi::RecvBuffer &) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}


template< typename ValueType >
void EdgeDoFAdditivePackInfo< ValueType >::communicateLocalCellToVertex(const Cell *, Vertex *) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}

///@}

template class EdgeDoFAdditivePackInfo< double >;
template class EdgeDoFAdditivePackInfo< float >;
template class EdgeDoFAdditivePackInfo< int >;

} // namespace hhg
