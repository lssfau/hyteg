#pragma once

#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

namespace hhg {

template< typename ValueType >
class BubblePackInfo : public communication::DoFSpacePackInfo< ValueType > {

public:
  BubblePackInfo(uint_t level,
                   PrimitiveDataID<VertexBubbleFunctionMemory< ValueType >, Vertex> dataIDVertex,
                   PrimitiveDataID<EdgeBubbleFunctionMemory< ValueType >, Edge> dataIDEdge,
                   PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> dataIDFace,
                   std::weak_ptr<PrimitiveStorage> storage)
    : communication::DoFSpacePackInfo< ValueType >(level, dataIDVertex, dataIDEdge, dataIDFace, storage){

  }
  virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  virtual void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const override;

  virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const override;

  virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const override;

  virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const override;

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
void BubblePackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

template< typename ValueType >
void BubblePackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

template< typename ValueType >
void BubblePackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
}

///@}
/// @name Edge to Vertex
///@{

template< typename ValueType >
void BubblePackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const {
  typedef stencilDirection sD;
  using namespace BubbleEdge;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  if(vertexIdOnEdge == 0){
    buffer << edgeData[edge_index(level_,0u,sD::CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      buffer << edgeData[edge_index(level_, 0u, sD::CELL_GRAY_NE)];
    }
  } else if(vertexIdOnEdge == 1){
    const uint_t lastPos = levelinfo::num_microvertices_per_edge(level_) -1;
    buffer << edgeData[edge_index(level_,lastPos,sD::CELL_GRAY_SW)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      buffer << edgeData[edge_index(level_, lastPos, sD::CELL_GRAY_NW)];
    }
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

template< typename ValueType >
void BubblePackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  for(PrimitiveID faceId : storage_.lock()->getEdge(sender)->neighborFaces())
  {
    const uint_t faceIdOnVertex = receiver->face_index(faceId);
    buffer >> vertexData[faceIdOnVertex];
  }
}

template< typename ValueType >
void BubblePackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const {
  using namespace BubbleEdge;
  typedef stencilDirection sD;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *vertexData = receiver->getData(dataIDVertex_)->getPointer( level_ );
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());

  if(vertexIdOnEdge == 0){
    const uint_t southFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[0]);
    vertexData[southFaceIdOnVertex] = edgeData[edge_index(level_,0u,sD::CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      const uint_t northFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[1]);
      vertexData[northFaceIdOnVertex] = edgeData[edge_index(level_, 0u, sD::CELL_GRAY_NE)];
    }
  } else if(vertexIdOnEdge == 1){
    const uint_t lastPos = levelinfo::num_microvertices_per_edge(level_) -1;
    const uint_t southFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[0]);
    vertexData[southFaceIdOnVertex] = edgeData[edge_index(level_,lastPos,sD::CELL_GRAY_SW)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      const uint_t northFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[1]);
      vertexData[northFaceIdOnVertex] = edgeData[edge_index(level_, lastPos, sD::CELL_GRAY_NW)];
    }
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver->getID().getID() << " is not in Edge: " << sender);
  }

}


///@}
/// @name Edge to Face
///@{

template< typename ValueType >
void BubblePackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(buffer);
}

template< typename ValueType >
void BubblePackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

template< typename ValueType >
void BubblePackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
}

///@}
/// @name Face to Edge
///@{

template< typename ValueType >
void BubblePackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const {
  using namespace hhg::BubbleFace;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_); it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_); it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
}

template< typename ValueType >
void BubblePackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const {
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t vPerEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t edgeIdOnFace = receiver->face_index(sender);
  stencilDirection dirCellGray;
  stencilDirection dirCellBlue;
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dirCellGray = stencilDirection::CELL_GRAY_SE;
    dirCellBlue = stencilDirection::CELL_BLUE_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dirCellGray = stencilDirection::CELL_GRAY_NE;
    dirCellBlue = stencilDirection::CELL_BLUE_NW;
  }
  //unpack Gray Cell
  for (uint_t i = 0; i < vPerEdge - 1; ++i)
  {
    buffer >> edgeData[BubbleEdge::edge_index(level_,i,dirCellGray)];
  }
  //unpack Blue Cell
  for (uint_t i = 1; i < vPerEdge - 1; ++i)
  {
    buffer >> edgeData[BubbleEdge::edge_index(level_,i,dirCellBlue)];
  }
}

template< typename ValueType >
void BubblePackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const {
  using namespace hhg::BubbleFace;
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIdOnFace = receiver->face_index(sender->getID());
  stencilDirection dirCellGray;
  stencilDirection dirCellBlue;
  uint_t edgeIndexOnFace = sender->edge_index(receiver->getID());
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dirCellGray = stencilDirection::CELL_GRAY_SE;
    dirCellBlue = stencilDirection::CELL_BLUE_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dirCellGray = stencilDirection::CELL_GRAY_NE;
    dirCellBlue = stencilDirection::CELL_BLUE_NW;
  }
  uint_t pos = 0;
  //copy Gray Cell
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_); it != indexIterator(); ++it)
  {
    edgeData[BubbleEdge::edge_index(level_,pos,dirCellGray)] = faceData[*it];
    pos++;
  }
  pos = 1;

  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_); it != indexIterator(); ++it)
  {
    edgeData[BubbleEdge::edge_index(level_,pos,dirCellBlue)] = faceData[*it];
    pos++;
  }
}

///@}

} //namespace hhg
