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
#include "N1E1PackInfo.hpp"

#include "hyteg/Algorithms.hpp"
#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {
namespace n1e1 {

template < typename ValueType >
N1E1PackInfo< ValueType >::N1E1PackInfo( uint_t                                                 level,
                                         PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                         PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                         PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                         PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                                         std::weak_ptr< PrimitiveStorage >                      storage )
: communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
, dofPackInfo_( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
{}

// TODO other functions
template < typename ValueType >
void N1E1PackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
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

   ValueType signX  = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;
   ValueType signY  = ( iterationVertex2 < iterationVertex0 ) ? -1 : 1;
   ValueType signXY = ( iterationVertex2 < iterationVertex1 ) ? -1 : 1;

   for ( const auto& cellIterator :
         edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
   {
      ValueType tmpX, tmpY, tmpXY;
      buffer >> tmpX;
      if ( globalDefines::useGeneratedKernels && level_ >= 1 )
      {
         buffer >> tmpXY;
         buffer >> tmpY;
      }
      else
      {
         buffer >> tmpY;
         buffer >> tmpXY;
      }

      cellData[edgedof::macrocell::index( level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationX )] =
          signX * tmpX;
      cellData[edgedof::macrocell::index( level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationY )] =
          signY * tmpY;
      cellData[edgedof::macrocell::index( level_, cellIterator.x(), cellIterator.y(), cellIterator.z(), dstEdgeOrientationXY )] =
          signXY * tmpXY;
   }
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell (unpack)" );
}

template < typename ValueType >
void N1E1PackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
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

   ValueType signX  = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;
   ValueType signY  = ( iterationVertex2 < iterationVertex0 ) ? -1 : 1;
   ValueType signXY = ( iterationVertex2 < iterationVertex1 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( const auto& faceIdx : edgedof::macroface::Iterator( level_ ) )
   {
      auto cellIdx = *cellIterator;

      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationX )] =
          signX * faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::X )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationY )] =
          signY * faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::Y )];
      cellData[edgedof::macrocell::index( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), dstEdgeOrientationXY )] =
          signXY * faceData[edgedof::macroface::index( level_, faceIdx.x(), faceIdx.y(), edgedof::EdgeDoFOrientation::XY )];

      cellIterator++;
   }

   WALBERLA_ASSERT( cellIterator == cellIterator.end() );
   this->storage_.lock()->getTimingTree()->stop( "EdgeDoF - Face to Cell" );
}

template < typename ValueType >
void N1E1PackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                    const PrimitiveID&         sender,
                                                    walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Cell -> Edge only meaningful in 3D." );

   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = receiver->getLocalEdgeID( sender );
   const uint_t iterationVertex0 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( { iterationVertex0, iterationVertex1 } ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType sign = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   for ( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level_ ); i++ )
   {
      ValueType tmp;
      buffer >> tmp;
      cellData[edgedof::macrocell::index( level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )] =
          sign * tmp;
      cellIterator++;
   }
}

template < typename ValueType >
void N1E1PackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Cell -> Edge only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver->getID() ) );

   ValueType* cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

   const uint_t localEdgeID      = receiver->getLocalEdgeID( sender->getID() );
   const uint_t iterationVertex0 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t iterationVertex1 = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   const uint_t iterationVertex2 =
       algorithms::getMissingIntegersAscending< 2, 4 >( { iterationVertex0, iterationVertex1 } ).at( 2 );

   auto srcEdgeOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
       edgedof::EdgeDoFOrientation::X, iterationVertex0, iterationVertex1, iterationVertex2 );

   ValueType sign = ( iterationVertex1 < iterationVertex0 ) ? -1 : 1;

   auto cellIterator = edgedof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

   const ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );

   for ( const auto& it : edgedof::macroedge::Iterator( level_ ) )
   {
      cellData[edgedof::macrocell::index( level_, cellIterator->x(), cellIterator->y(), cellIterator->z(), srcEdgeOrientation )] =
          sign * edgeData[edgedof::macroedge::index( level_, it.x() )];
      cellIterator++;
   }
}

template class N1E1PackInfo< double >;
template class N1E1PackInfo< int >;
template class N1E1PackInfo< long >;
template class N1E1PackInfo< long long >;

} // namespace n1e1
} // namespace hyteg
