#pragma once

#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

namespace hhg {

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
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 0 ) )
  {
    buffer >> faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  uint_t idx = 0;
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 0 ) )
  {
    faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ] = edgeData[idx];
    idx++;
  }
}

///@}
/// @name Face to Edge
///@{

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 1 ) )
  {
    buffer << faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t edgeIdOnFace = receiver->face_index(sender);
  stencilDirection dir;
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dir = stencilDirection::VERTEX_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dir = stencilDirection::VERTEX_N;
  }
  for (uint_t i = 0; i < v_perEdge -1; ++i)
  {
    buffer >> edgeData[ vertexdof::macroedge::indexFromVertex( level_, i, dir ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const
{
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  stencilDirection dir;
  uint_t edgeIdOnFace = sender->edge_index(receiver->getID());
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dir = stencilDirection::VERTEX_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dir = stencilDirection::VERTEX_N;
  }
  uint_t idx = 0;
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIdOnFace, sender->edge_orientation[edgeIdOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 1 ) )
  {
    edgeData[ vertexdof::macroedge::indexFromVertex( level_, idx, dir ) ] = faceData[ vertexdof::macroface::indexFromVertex( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
    idx++;
  }
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

  for ( const auto & it : vertexdof::macrocell::BorderIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 ) )
  {
    buffer >> cellData[ vertexdof::macrocell::indexFromVertex( level_, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalFaceToCell(const Face *sender, Cell *receiver) const
{
  const ValueType * faceData = sender->getData( dataIDFace_ )->getPointer( level_ );
        ValueType * cellData = receiver->getData( dataIDCell_ )->getPointer( level_ );

  const uint_t localFaceID = receiver->getLocalFaceID( sender->getID() );
  const uint_t iterationVertex0 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = receiver->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  auto cellIterator = vertexdof::macrocell::BorderIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2 );

  for ( const auto & faceIdx : vertexdof::macroface::Iterator( level_ ) )
  {
    auto cellIdx = *cellIterator;

    cellData[ vertexdof::macrocell::indexFromVertex( level_, cellIdx.x(), cellIdx.y(), cellIdx.z(), stencilDirection::VERTEX_C ) ] =
        faceData[ vertexdof::macroface::indexFromVertex( level_, faceIdx.x(), faceIdx.y(), stencilDirection::VERTEX_C ) ];

    cellIterator++;
  }

  WALBERLA_ASSERT( cellIterator == cellIterator.end() );
}

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::packCellForFace(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  const ValueType * cellData = sender->getData( dataIDCell_ )->getPointer( level_ );
  const uint_t localFaceID = sender->getLocalFaceID( receiver );
  const uint_t iterationVertex0 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 );
  const uint_t iterationVertex1 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 );
  const uint_t iterationVertex2 = sender->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 );

  for ( const auto & it : vertexdof::macrocell::BorderIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 ) )
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

template< typename ValueType >
void VertexDoFPackInfo< ValueType >::communicateLocalCellToFace(const Cell *sender, Face *receiver) const
{
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

  auto cellIterator = vertexdof::macrocell::BorderIterator( level_, iterationVertex0, iterationVertex1, iterationVertex2, 1 );

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
}



} //namespace hhg
