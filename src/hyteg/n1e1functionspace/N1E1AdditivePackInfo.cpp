/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "N1E1AdditivePackInfo.hpp"

#include "hyteg/Algorithms.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {
namespace n1e1 {

template < typename ValueType >
N1E1AdditivePackInfo< ValueType >::N1E1AdditivePackInfo( uint_t                                                 level,
                                                         PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                                         PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                                         PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                                         PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                                                         std::weak_ptr< PrimitiveStorage >                      storage )
: communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
, dofPackInfo_( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
{}

/// @name Vertex to Edge
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packVertexForEdge( const Vertex*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackEdgeFromVertex( Edge*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex*, Edge* ) const
{
   WALBERLA_ABORT( "Additive communication Vertex -> Edge not supported." );
}

///@}
/// @name Edge to Vertex
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packEdgeForVertex( const Edge*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackVertexFromEdge( Vertex*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge*, Vertex* ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Vertex not supported." );
}

///@}
/// @name Edge to Face
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packEdgeForFace( const Edge*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackFaceFromEdge( Face*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalEdgeToFace( const Edge*, Face* ) const
{
   WALBERLA_ABORT( "Additive communication Edge -> Face not supported." );
}

///@}
/// @name Face to Edge
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packFaceForEdge( const Face*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );
   WALBERLA_ABORT( "Not implemented." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackEdgeFromFace( Edge*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );
   WALBERLA_ABORT( "Not implemented." );
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalFaceToEdge( const Face*, Edge* ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Edge only meaningful in 2D." );
   WALBERLA_ABORT( "Not implemented." );
}

///@}
/// @name Face to Vertex
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packFaceForVertex( const Face*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackVertexFromFace( Vertex*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalFaceToVertex( const Face*, Vertex* ) const
{
   WALBERLA_CHECK( !this->storage_.lock()->hasGlobalCells(), "Additive communication Face -> Vertex only meaningful in 2D." );
   // nothing to be done
}

///@}
/// @name Cell to Face
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packCellForFace( const Cell*                sender,
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

   ValueType signX  = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;
   ValueType signY  = ( iterationVertex2 < iterationVertex0 ) ? -1 : 1;
   ValueType signXY = ( iterationVertex2 < iterationVertex1 ) ? -1 : 1;

   for ( const auto& cellIterator :
         edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      buffer << signX * cellData[edgedof::macrocell::index(
                            level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationX )];
      buffer << signY * cellData[edgedof::macrocell::index(
                            level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationY )];
      buffer << signXY * cellData[edgedof::macrocell::index(
                             level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), srcEdgeOrientationXY )];
   }
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Face only meaningful in 3D." );

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

   ValueType signX  = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;
   ValueType signY  = ( iterationVertex2 < iterationVertex0 ) ? -1 : 1;
   ValueType signXY = ( iterationVertex2 < iterationVertex1 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )] +=
          signX * cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationX )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )] +=
          signY * cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationY )];
      faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )] +=
          signXY * cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), srcEdgeOrientationXY )];
      cellIterator++;
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
}

///@}
/// @name Cell to Edge
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packCellForEdge( const Cell*                sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = sender->getLocalEdgeID( receiver );
   const uint_t iterationVertex0 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( { iterationVertex0, iterationVertex1 } ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType sign = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); i++ )
   {
      buffer << sign * cellData[edgedof::macrocell::index(
                           level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )];
      cellIterator++;
   }
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalCellToEdge( const Cell* sender, Edge* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

   ValueType* cellData = sender->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = sender->getLocalEdgeID( receiver->getID() );
   const uint_t iterationVertex0 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = sender->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( { iterationVertex0, iterationVertex1 } ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType sign = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType* edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );

   for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
   {
      edgeData[edgedof::macroedge::index( level_, it.x() )] +=
          sign * cellData[edgedof::macrocell::index(
                     level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )];
      cellIterator++;
   }
}

///@}
/// @name Cell to Vertex
///@{

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::packCellForVertex( const Cell*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::unpackVertexFromCell( Vertex*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}

template < typename ValueType >
void N1E1AdditivePackInfo< ValueType >::communicateLocalCellToVertex( const Cell*, Vertex* ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Additive communication Cell -> Vertex only meaningful in 3D." );
   // nothing to be done
}

///@}

template class N1E1AdditivePackInfo< double >;
template class N1E1AdditivePackInfo< int >;
template class N1E1AdditivePackInfo< long >;
template class N1E1AdditivePackInfo< long long >;

} // namespace n1e1
} // namespace hyteg
