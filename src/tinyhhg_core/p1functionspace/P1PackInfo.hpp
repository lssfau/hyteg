#pragma once

#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg {

namespace {
SPECIALIZE( uint_t, vertexdof::macroedge::indexFromVertex, indexFromVertexOnMacroEdge );
SPECIALIZE( uint_t, vertexdof::macroface::indexFromVertex, indexFromVertexOnMacroFace );
}

template< typename ValueType >
class P1PackInfo : public communication::DoFSpacePackInfo<ValueType> {

public:
  P1PackInfo(uint_t level,
             PrimitiveDataID<FunctionMemory< ValueType >, Vertex> dataIDVertex,
             PrimitiveDataID<FunctionMemory< ValueType >, Edge> dataIDEdge,
             PrimitiveDataID<FunctionMemory< ValueType >, Face> dataIDFace,
             std::weak_ptr<PrimitiveStorage> storage)
    : communication::DoFSpacePackInfo< ValueType >(level,dataIDVertex,dataIDEdge,dataIDFace,storage)
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

private:
  using communication::DoFSpacePackInfo< ValueType >::level_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
  using communication::DoFSpacePackInfo< ValueType >::storage_;
};

/// @name Vertex to Edge
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const {
  WALBERLA_UNUSED(receiver);
  ValueType *vertexData = sender->getData(dataIDVertex_)->getPointer( level_ );
  buffer << vertexData[0];
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
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
  buffer >> edgeData[indexFromVertexOnMacroEdge( level_, pos, stencilDirection::VERTEX_C ) ];
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const
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
  edgeData[indexFromVertexOnMacroEdge( level_, pos, stencilDirection::VERTEX_C ) ] = vertexData[0];
}

///@}
/// @name Edge to Vertex
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    buffer << edgeData[indexFromVertexOnMacroEdge( level_, 1u, stencilDirection::VERTEX_C ) ];
  } else if(vertexIdOnEdge == 1){
    buffer << edgeData[indexFromVertexOnMacroEdge( level_, levelinfo::num_microvertices_per_edge(level_)-2 ,stencilDirection::VERTEX_C ) ];
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
{
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t edgeIdOnVertex = receiver->edge_index(sender);
  buffer >> vertexData[edgeIdOnVertex + 1];
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());
  uint_t edgeIdOnVertex = receiver->edge_index(sender->getID());
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    const uint_t idx = indexFromVertexOnMacroEdge( level_, 1u, stencilDirection::VERTEX_C );
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else if(vertexIdOnEdge == 1){
    const uint_t idx = indexFromVertexOnMacroEdge( level_, levelinfo::num_microvertices_per_edge(level_)-2, stencilDirection::VERTEX_C );
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << receiver << " is not contained in Edge: " << sender);
  }
}

///@}
/// @name Edge to Face
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);

  for (uint_t i = 0; i < v_perEdge; ++i) {
    buffer << edgeData[ indexFromVertexOnMacroEdge( level_, i, stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  using namespace hhg::P1Face;
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 0 ) )
  {
    buffer >> faceData[ indexFromVertexOnMacroFace( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const
{
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  uint_t idx = 0;
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 0 ) )
  {
    faceData[ indexFromVertexOnMacroFace( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ] = edgeData[idx];
    idx++;
  }
}

///@}
/// @name Face to Edge
///@{

template< typename ValueType >
void P1PackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const
{
  using namespace hhg::P1Face;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  indexing::FaceBorderDirection faceBorderDirection = indexing::getFaceBorderDirection( edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace] );
  for( const auto & it : vertexdof::macroface::BorderIterator( level_, faceBorderDirection, 1 ) )
  {
    buffer << faceData[ indexFromVertexOnMacroFace( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const
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
    buffer >> edgeData[ indexFromVertexOnMacroEdge( level_, i, dir ) ];
  }
}

template< typename ValueType >
void P1PackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const
{
  using namespace hhg::P1Face;
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
    edgeData[ indexFromVertexOnMacroEdge( level_, idx, dir ) ] = faceData[ indexFromVertexOnMacroFace( level_, it.col(), it.row(), stencilDirection::VERTEX_C ) ];
    idx++;
  }
}


} //namespace hhg
