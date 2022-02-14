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
#include "hyteg/primitives/all.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/generatedKernels/communicate_directly_vertexdof_cell_to_face.hpp"
#include "hyteg/p1functionspace/generatedKernels/communicate_directly_vertexdof_face_to_cell.hpp"

namespace hyteg {

template< typename ValueType >
class VertexDoFPackInfo : public communication::DoFSpacePackInfo<ValueType>
{
public:

  VertexDoFPackInfo( uint_t level,
                     PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                     PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                     PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                     std::weak_ptr< PrimitiveStorage>                       storage )
    : communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, storage )
  {}

  VertexDoFPackInfo(uint_t level,
             PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
             PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
             PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
             PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
             std::weak_ptr< PrimitiveStorage >                      storage )
    : communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage )
  {}

  void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const override;

  void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const override;

  void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const override;

  void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const override;

  void packFaceForCell(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackCellFromFace(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalFaceToCell(const Face *sender, Cell *receiver) const override;

  void packCellForFace(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackFaceFromCell(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalCellToFace(const Cell *sender, Face *receiver) const override;

  void packVertexForCell(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackCellFromVertex(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalVertexToCell(const Vertex *sender, Cell *receiver) const override;
  
  void packEdgeForCell(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackCellFromEdge(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalEdgeToCell(const Edge *sender, Cell *receiver) const override;



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

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const {
  WALBERLA_UNUSED(receiver);
  ValueType *vertexData = sender->getData(dataIDVertex_)->getPointer( level_ );
  buffer << vertexData[0];
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType* edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  //position in edge memory
  uint_t pos;
  if(receiver->vertex_index(sender) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << sender.getID() << " is not in Edge: " << receiver)
  }
  buffer >> edgeData[vertexdof::macroedge::indexFromVertex( level_, pos, stencilDirection::VERTEX_C ) ];
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const
{
  ValueType *vertexData = sender->getData(dataIDVertex_)->getPointer( level_ );
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t pos;
  if(receiver->vertex_index(sender->getID()) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender->getID()) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << sender << " is not in Edge: " << receiver)
  }
  edgeData[vertexdof::macroedge::indexFromVertex( level_, pos, stencilDirection::VERTEX_C ) ] = vertexData[0];
}

///@}
/// @name Edge to Vertex
///@{

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    buffer << edgeData[vertexdof::macroedge::indexFromVertex( level_, 1u, stencilDirection::VERTEX_C ) ];
  } else if(vertexIdOnEdge == 1){
    buffer << edgeData[vertexdof::macroedge::indexFromVertex( level_, levelinfo::num_microvertices_per_edge(level_)-2 ,stencilDirection::VERTEX_C ) ];
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t edgeIdOnVertex = receiver->edge_index(sender);
  buffer >> vertexData[edgeIdOnVertex + 1];
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());
  uint_t edgeIdOnVertex = receiver->edge_index(sender->getID());
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    const uint_t idx = vertexdof::macroedge::indexFromVertex( level_, 1u, stencilDirection::VERTEX_C );
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else if(vertexIdOnEdge == 1){
    const uint_t idx = vertexdof::macroedge::indexFromVertex( level_, levelinfo::num_microvertices_per_edge(level_)-2, stencilDirection::VERTEX_C );
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << receiver << " is not contained in Edge: " << sender);
  }
}

///@}
/// @name Edge to Face
///@{

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);

  for (uint_t i = 0; i < v_perEdge; ++i) {
    buffer << edgeData[ vertexdof::macroedge::indexFromVertex( level_, i, stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  indexing::FaceBoundaryDirection faceBorderDirection =
      indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0 ) )
  {
    buffer >> faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Edge to Face" );
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  uint_t idx = 0;
  indexing::FaceBoundaryDirection faceBorderDirection =
      indexing::getFaceBoundaryDirection( edgeIndexOnFace, receiver->getEdgeOrientation()[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 0 ) )
  {
    faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ] = edgeData[idx];
    idx++;
  }
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Edge to Face" );
}

///@}
/// @name Face to Edge
///@{

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  indexing::FaceBoundaryDirection faceBorderDirection =
      indexing::getFaceBoundaryDirection( edgeIndexOnFace, sender->getEdgeOrientation()[edgeIndexOnFace] );

  for( const auto & it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 1 ) )
  {
    buffer << faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }

  // To pack DoFs on face ghost-layers, we use an iterator with a width reduced by 1
  // Top ghost-layer is sent if there are any neighboring cells
  if ( sender->getNumNeighborCells() > 0 )
  {
    for ( const auto & it : indexing::FaceBoundaryIterator( levelinfo::num_microvertices_per_edge( level_ ) - 1, faceBorderDirection, 1 ))
    {
      buffer << faceData[vertexdof::macroface::index( level_, it.col(), it.row(), 0 )];
    }

    // Bottom ghost-layer is only sent if there is a second neighboring cell
    if ( sender->getNumNeighborCells() == 2 )
    {
      for ( const auto & it : indexing::FaceBoundaryIterator( levelinfo::num_microvertices_per_edge( level_ ) - 1, faceBorderDirection, 1 ))
      {
        buffer << faceData[vertexdof::macroface::index( level_, it.col(), it.row(), 1 )];
      }
    }
  }

}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  const uint_t faceIDOnEdge = receiver->face_index( sender );
  const auto storage = storage_.lock();
  WALBERLA_CHECK_NOT_NULLPTR( storage.get() );
  WALBERLA_CHECK( storage->faceExistsLocally( sender ) || storage->faceExistsInNeighborhood( sender ) );
  const auto senderFace = storage->getFace( sender );

  for (uint_t i = 0; i < vertexdof::macroedge::neighborFaceGhostLayerSize( level_ ); ++i)
  {
    buffer >> edgeData[ vertexdof::macroedge::indexOnNeighborFace( level_, i, faceIDOnEdge ) ];
  }

  // Unpacking the DoFs from the face ghost-layers (located in the interior of a macro-cell) now.
  if ( senderFace->getNumNeighborCells() > 0 )
  {
    const auto topCellPrimitiveID = senderFace->getCellID0();
    const auto localTopCellIDOnEdge = receiver->cell_index( topCellPrimitiveID );
    for (uint_t i = 0; i < vertexdof::macroedge::neighborCellGhostLayerSize( level_ ); ++i)
    {
      buffer >> edgeData[ vertexdof::macroedge::indexOnNeighborCell( level_, i, localTopCellIDOnEdge, receiver->getNumNeighborFaces() ) ];
    }

    if ( senderFace->getNumNeighborCells() == 2 )
    {
      const auto bottomCellPrimitiveID = senderFace->getCellID1();
      const auto localBottomCellIDOnEdge = receiver->cell_index( bottomCellPrimitiveID );
      for (uint_t i = 0; i < vertexdof::macroedge::neighborCellGhostLayerSize( level_ ); ++i)
      {
        buffer >> edgeData[ vertexdof::macroedge::indexOnNeighborCell( level_, i, localBottomCellIDOnEdge, receiver->getNumNeighborFaces() ) ];
      }
    }
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Face to Edge" );
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  uint_t edgeIdOnFace = sender->edge_index(receiver->getID());
  indexing::FaceBoundaryDirection faceBorderDirection =
      indexing::getFaceBoundaryDirection( edgeIdOnFace, sender->getEdgeOrientation()[edgeIdOnFace] );
  uint_t idx = 0;
  for( const auto & it : vertexdof::macroface::BoundaryIterator( level_, faceBorderDirection, 1 ) )
  {
    edgeData[ vertexdof::macroedge::indexOnNeighborFace( level_, idx, faceIdOnEdge ) ] = faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
    idx++;
  }

  // To pack DoFs on face ghost-layers, we use an iterator with a width reduced by 1
  // Top ghost-layer is sent if there are any neighboring cells
  if ( sender->getNumNeighborCells() > 0 )
  {
    const auto topCellPrimitiveID = sender->getCellID0();
    const auto localTopCellIDOnEdge = receiver->cell_index( topCellPrimitiveID );
    idx = 0;
    for ( const auto & it : indexing::FaceBoundaryIterator( levelinfo::num_microvertices_per_edge( level_ ) - 1, faceBorderDirection, 1 ))
    {
      edgeData[ vertexdof::macroedge::indexOnNeighborCell( level_, idx, localTopCellIDOnEdge, receiver->getNumNeighborFaces() ) ] = faceData[vertexdof::macroface::index( level_, it.col(), it.row(), 0 )];
      idx++;
    }
    // Bottom ghost-layer is only sent if there is a second neighboring cell
    if ( sender->getNumNeighborCells() == 2 )
    {
      const auto bottomCellPrimitiveID = sender->getCellID1();
      const auto localBottomCellIDOnEdge = receiver->cell_index( bottomCellPrimitiveID );
      idx = 0;
      for ( const auto & it : indexing::FaceBoundaryIterator( levelinfo::num_microvertices_per_edge( level_ ) - 1, faceBorderDirection, 1 ))
      {
        edgeData[ vertexdof::macroedge::indexOnNeighborCell( level_, idx, localBottomCellIDOnEdge, receiver->getNumNeighborFaces() ) ] = faceData[vertexdof::macroface::index( level_, it.col(), it.row(), 1 )];
        idx++;
      }
    }
  }
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Face to Edge" );
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packFaceForCell(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  const ValueType * faceData = sender->getData( dataIDFace_ )->getPointer( level_ );

  // only inner points
  for ( const auto & it : vertexdof::macroface::Iterator( level_ ) )
  {
    buffer << faceData[ vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackCellFromFace(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType * cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );
  const uint_t localFaceID = receiver->getLocalFaceID( sender );
  const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
  {
    buffer >> cellData[ vertexdof::macrocell::indexFromVertex( level_, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ];
  }
}

template<>
inline void VertexDoFPackInfo< real_t >::communicateLocalFaceToCell(const Face *sender, Cell *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Face to Cell" );
  const real_t * faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
  real_t * cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

  const uint_t localFaceID = receiver->getLocalFaceID( sender->getID());
  const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  if ( globalDefines::useGeneratedKernels )
  {
     vertexdof::comm::generated::communicate_directly_vertexdof_face_to_cell( cellData,
                                                                              faceData,
                                                                              static_cast< int32_t >( level_ ),
                                                                              static_cast< int64_t >( iterationVertex0 ),
                                                                              static_cast< int64_t >( iterationVertex1 ),
                                                                              static_cast< int64_t >( iterationVertex2 ) );
  }
  else
  {
    auto cellIterator = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

    for ( const auto & faceIdx : vertexdof::macroface::Iterator( level_ ))
    {
      auto cellIdx = *cellIterator;

      cellData[vertexdof::macrocell::indexFromVertex( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C )] =
      faceData[vertexdof::macroface::indexFromVertex( level_, faceIdx.x(), faceIdx.y(), stencilDirection::VERTEX_C )];

      cellIterator++;
    }

    WALBERLA_ASSERT( cellIterator == cellIterator.end());
  }
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Face to Cell" );
}

template< typename ValueType >
inline void VertexDoFPackInfo< ValueType >::communicateLocalFaceToCell(const Face *sender, Cell *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Face to Cell" );
  const ValueType * faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
        ValueType * cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

  const uint_t localFaceID = receiver->getLocalFaceID( sender->getID() );
  const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  auto cellIterator = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

  for ( const auto & faceIdx : vertexdof::macroface::Iterator( level_ ) )
  {
    auto cellIdx = *cellIterator;

    cellData[ vertexdof::macrocell::indexFromVertex( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C ) ] =
        faceData[ vertexdof::macroface::indexFromVertex( level_, faceIdx.x(), faceIdx.y(), stencilDirection::VERTEX_C ) ];

    cellIterator++;
  }

  WALBERLA_ASSERT( cellIterator == cellIterator.end() );
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Face to Cell" );
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packCellForFace(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  const ValueType * cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
  const uint_t localFaceID = sender->getLocalFaceID( receiver );
  const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 ) )
  {
    buffer << cellData[ vertexdof::macrocell::indexFromVertex( level_, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ];
  }
}


template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackFaceFromCell(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType * faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

  stencilDirection neighborDirection;

  WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
  WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender ) );

  if ( receiver->cell_index( sender ) == 0 )
  {
    neighborDirection = stencilDirection::VERTEX_TC;
  }
  else
  {
    WALBERLA_ASSERT_EQUAL( receiver->cell_index( sender ), 1 );
    neighborDirection = stencilDirection::VERTEX_BC;
  }

  for ( const auto & it : vertexdof::macroface::Iterator( level_ ) )
  {
    if ( it.x() + it.y() < levelinfo::num_microvertices_per_edge( level_ ) - 1 )
    {
      buffer >> faceData[ vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), neighborDirection ) ];
    }
  }
}


template<>
inline void VertexDoFPackInfo< real_t >::communicateLocalCellToFace(const Cell *sender, Face *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Cell to Face" );
  const real_t * cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
  const uint_t localFaceID = sender->getLocalFaceID( receiver->getID() );
  const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  real_t * faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

  WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
  WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

  if ( globalDefines::useGeneratedKernels )
  {
    const auto faceLocalCellID = receiver->cell_index( sender->getID() );
    const auto offsetToGhostLayer =
        faceLocalCellID == 0 ?
            levelinfo::num_microvertices_per_face( level_ ) :
            levelinfo::num_microvertices_per_face( level_ ) +
                levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level_ ) - 1 );
    vertexdof::comm::generated::communicate_directly_vertexdof_cell_to_face( cellData,
                                                                             &faceData[offsetToGhostLayer],
                                                                             static_cast< int32_t >( level_ ),
                                                                             static_cast< int64_t >( iterationVertex0 ),
                                                                             static_cast< int64_t >( iterationVertex1 ),
                                                                             static_cast< int64_t >( iterationVertex2 ) );
  }
  else
  {
    stencilDirection neighborDirection;

    if ( receiver->cell_index( sender->getID()) == 0 )
    {
      neighborDirection = stencilDirection::VERTEX_TC;
    } else
    {
      WALBERLA_ASSERT_EQUAL( receiver->cell_index( sender->getID()), 1 );
      neighborDirection = stencilDirection::VERTEX_BC;
    }

    auto cellIterator = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );

    for ( const auto & it : vertexdof::macroface::Iterator( level_ ))
    {
      if ( it.x() + it.y() < levelinfo::num_microvertices_per_edge( level_ ) - 1 )
      {
        auto cellIdx = *cellIterator;
        faceData[vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), neighborDirection )] =
        cellData[vertexdof::macrocell::indexFromVertex( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C )];
        cellIterator++;
      }
    }

    WALBERLA_ASSERT( cellIterator == cellIterator.end());
  }
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Cell to Face" );
}


template< typename ValueType >
inline void VertexDoFPackInfo< ValueType >::communicateLocalCellToFace(const Cell *sender, Face *receiver) const
{
  this->storage_.lock()->getTimingTree()->start( "VertexDoF - Cell to Face" );
  const ValueType * cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
  const uint_t localFaceID = sender->getLocalFaceID( receiver->getID() );
  const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  ValueType * faceData = receiver->getData( dataIDFace_ )->getPointer( level_ );

  stencilDirection neighborDirection;

  WALBERLA_ASSERT_GREATER( receiver->getNumNeighborCells(), 0 );
  WALBERLA_ASSERT( receiver->neighborPrimitiveExists( sender->getID() ) );

  if ( receiver->cell_index( sender->getID() ) == 0 )
  {
    neighborDirection = stencilDirection::VERTEX_TC;
  }
  else
  {
    WALBERLA_ASSERT_EQUAL( receiver->cell_index( sender->getID() ), 1 );
    neighborDirection = stencilDirection::VERTEX_BC;
  }

  auto cellIterator = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );

  for ( const auto & it : vertexdof::macroface::Iterator( level_ ) )
  {
    if ( it.x() + it.y() < levelinfo::num_microvertices_per_edge( level_ ) - 1 )
    {
      auto cellIdx = *cellIterator;
      faceData[ vertexdof::macroface::indexFromVertex( level_, it.x(), it.y(), neighborDirection ) ] =
          cellData[ vertexdof::macrocell::indexFromVertex( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C ) ];
      cellIterator++;
    }
  }

  WALBERLA_ASSERT( cellIterator == cellIterator.end() );
  this->storage_.lock()->getTimingTree()->stop( "VertexDoF - Cell to Face" );
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::packVertexForCell( const Vertex*              sender,
                                                        const PrimitiveID&         receiver,
                                                        walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Vertex -> Cell only meaningful in 3D." );

   const ValueType* vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver ) );

   buffer << vertexData[0];
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackCellFromVertex( Cell*                      receiver,
                                                           const PrimitiveID&         sender,
                                                           walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Vertex -> Cell only meaningful in 3D." );

   ValueType*      cellData      = receiver->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t    localVertexID = receiver->getLocalVertexID( sender );
   indexing::Index microVertexIndexInMacroCell( 0, 0, 0 );
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
   ValueType tmp;
   buffer >> tmp;
   cellData[vertexdof::macrocell::indexFromVertex( level_,
                                                   microVertexIndexInMacroCell.x(),
                                                   microVertexIndexInMacroCell.y(),
                                                   microVertexIndexInMacroCell.z(),
                                                   stencilDirection::VERTEX_C )] = tmp;
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Vertex -> Cell only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver->getID() ) );

   const ValueType* vertexData    = sender->getData( dataIDVertex_ )->getPointer( level_ );
   ValueType*       cellData      = receiver->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t     localVertexID = receiver->getLocalVertexID( sender->getID() );
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
   cellData[vertexdof::macrocell::indexFromVertex( level_,
                                                   microVertexIndexInMacroCell.x(),
                                                   microVertexIndexInMacroCell.y(),
                                                   microVertexIndexInMacroCell.z(),
                                                   stencilDirection::VERTEX_C )] = vertexData[0];
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::packEdgeForCell( const Edge*                sender,
                                                      const PrimitiveID&         receiver,
                                                      walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Edge -> Cell only meaningful in 3D." );

   const ValueType* edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver ) );

   for ( const auto& it : vertexdof::macroedge::Iterator( level_ ) )
   {
      buffer << edgeData[vertexdof::macroedge::indexFromVertex( level_, it.x(), stencilDirection::VERTEX_C )];
   }
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                         const PrimitiveID&         sender,
                                                         walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Edge -> Cell only meaningful in 3D." );

   ValueType*         cellData                  = receiver->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t       localEdgeID               = receiver->getLocalEdgeID( sender );
   const uint_t       iterationVertex0          = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t       iterationVertex1          = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   std::set< uint_t > possibleIterationVertices = {0, 1, 2, 3};
   possibleIterationVertices.erase( iterationVertex0 );
   possibleIterationVertices.erase( iterationVertex1 );
   const uint_t iterationVertex2 = *possibleIterationVertices.begin();

   const uint_t edgeSize = levelinfo::num_microvertices_per_edge( level_ );
   auto         it = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 );
   for ( uint_t i = 0; i < edgeSize; i++ )
   {
      ValueType tmp;
      buffer >> tmp;
      cellData[vertexdof::macrocell::indexFromVertex( level_, it->x(), it->y(), it->z(), stencilDirection::VERTEX_C )] = tmp;
      it++;
   }
}

template < typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   WALBERLA_CHECK( this->storage_.lock()->hasGlobalCells(), "Communication Edge -> Cell only meaningful in 3D." );

   WALBERLA_ASSERT_GREATER( sender->getNumNeighborCells(), 0 );
   WALBERLA_ASSERT( sender->neighborPrimitiveExists( receiver->getID() ) );

   const ValueType*   edgeData                  = sender->getData( dataIDEdge_ )->getPointer( level_ );
   ValueType*         cellData                  = receiver->getData( dataIDCell_ )->getPointer( level_ );
   const uint_t       localEdgeID               = receiver->getLocalEdgeID( sender->getID() );
   const uint_t       iterationVertex0          = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 0 );
   const uint_t       iterationVertex1          = receiver->getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeID ).at( 1 );
   std::set< uint_t > possibleIterationVertices = {0, 1, 2, 3};
   possibleIterationVertices.erase( iterationVertex0 );
   possibleIterationVertices.erase( iterationVertex1 );
   const uint_t iterationVertex2 = *possibleIterationVertices.begin();

   const uint_t edgeSize = levelinfo::num_microvertices_per_edge( level_ );
   auto         it = vertexdof::macrocell::BoundaryIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 0 );
   for ( uint_t i = 0; i < edgeSize; i++ )
   {
      cellData[vertexdof::macrocell::indexFromVertex( level_, it->x(), it->y(), it->z(), stencilDirection::VERTEX_C )] =
          edgeData[i];
      it++;
   }
}

} //namespace hyteg
