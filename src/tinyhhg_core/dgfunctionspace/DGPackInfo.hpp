#pragma once

#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "DGEdgeIndex.hpp"

namespace hhg{

using walberla::uint_t;

template< typename ValueType >
class DGPackInfo : public communication::DoFSpacePackInfo< ValueType > {

public:
  DGPackInfo(uint_t level,
                 PrimitiveDataID<FunctionMemory< ValueType >, Vertex> dataIDVertex,
                 PrimitiveDataID<FunctionMemory< ValueType >, Edge> dataIDEdge,
                 PrimitiveDataID<FunctionMemory< ValueType >, Face> dataIDFace,
                 std::weak_ptr<PrimitiveStorage> storage)
      : communication::DoFSpacePackInfo< ValueType >(level, dataIDVertex, dataIDEdge, dataIDFace, storage){

  }
  void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) override;

  void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) override;

  void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) override;

  void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) override;

  void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) override;

  void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) override;

private:
  using communication::DoFSpacePackInfo< ValueType >::level_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
  using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
  using communication::DoFSpacePackInfo< ValueType >::storage_;

};

template< typename ValueType >
void DGPackInfo< ValueType >::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
/// see DGMemory.hpp for a description of the Vertex Memory
  ValueType *vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );
  for(const PrimitiveID& faceID: storage_.lock()->getEdge(receiver)->neighborFaces()) {
    buffer << vertexData[ sender->face_index(faceID) * 2];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  typedef stencilDirection sD;
  ValueType *edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t pos = std::numeric_limits<uint_t>::max();
  if(receiver->vertex_index(sender) == 0) {
    pos = 0;
  } else if (receiver->vertex_index(sender) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 2;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << sender.getID() << " is not in Edge: " << receiver)
  }
  buffer >> edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_GRAY_SE)];
  if(receiver->getNumNeighborFaces() == 2){
    buffer >> edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_GRAY_NE)];
  }

}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) {
  typedef stencilDirection sD;
  ValueType *vertexData = sender->getData( dataIDVertex_ )->getPointer( level_ );
  ValueType *edgeData = receiver->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t pos = std::numeric_limits<uint_t>::max();
  if(receiver->vertex_index(sender->getID()) == 0) {
    pos = 0;
  } else if (receiver->vertex_index(sender->getID()) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 2;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << sender << " is not in Edge: " << receiver)
  }
  edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_GRAY_SE)] = vertexData[ sender->face_index(receiver->neighborFaces()[0]) * 2];
  if(receiver->getNumNeighborFaces() == 2){
    edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_GRAY_NE)] = vertexData[ sender->face_index(receiver->neighborFaces()[1]) * 2];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  ///the blue face DoF which are owned by the face need to communicated to the vertex
  typedef stencilDirection sD;
  ValueType *edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t pos = std::numeric_limits<uint_t>::max();
  if(sender->vertex_index(receiver) == 0) {
    pos = 1;
  } else if (sender->vertex_index(receiver) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 2;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender)
  }
  buffer << edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_BLUE_SE)];
  if(sender-> getNumNeighborFaces() == 2){
    buffer << edgeData[BubbleEdge::edge_index(level_,pos,sD::CELL_BLUE_NW)];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  ValueType *vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
  for(const PrimitiveID& faceID: storage_.lock()->getEdge(sender)->neighborFaces()) {
    buffer >> vertexData[ receiver->face_index(faceID) * 2 + 1];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {
  ValueType *edgeData = sender->getData( dataIDEdge_ )->getPointer( level_ );
  uint_t pos = std::numeric_limits<uint_t>::max();
  if(sender->vertex_index(receiver->getID()) == 0) {
    pos = 1;
  } else if (sender->vertex_index(receiver->getID()) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 2;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver << " is not in Edge: " << sender)
  }
  ValueType *vertexData = receiver->getData( dataIDVertex_ )->getPointer( level_ );
  vertexData[ receiver->face_index(sender->neighborFaces()[0]) * 2 + 1] =
    edgeData[BubbleEdge::edge_index(level_,pos,stencilDirection::CELL_BLUE_SE)];;
  if(sender->getNumNeighborFaces() == 2){
    vertexData[ receiver->face_index(sender->neighborFaces()[1]) * 2 + 1] =
      edgeData[BubbleEdge::edge_index(level_,pos,stencilDirection::CELL_BLUE_NW)];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  uint_t vPerEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t faceIdOnEdge = sender->face_index(receiver);
  stencilDirection dirCellGray;
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dirCellGray = stencilDirection::CELL_GRAY_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dirCellGray = stencilDirection::CELL_GRAY_NE;
  }
  for (uint_t i = 0; i < vPerEdge - 1; ++i)
  {
    buffer << edgeData[BubbleEdge::edge_index(level_,i,dirCellGray)];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  using namespace BubbleFace;
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  for(auto it = indexIterator(edgeIndexOnFace,
                              receiver->edge_orientation[edgeIndexOnFace],
                              CELL_GRAY, level_);
      it != indexIterator();
      ++it){
    buffer >> faceData[*it];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {
  using namespace BubbleFace;
  ValueType *edgeData = sender->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = receiver->getData(dataIDFace_)->getPointer( level_ );
  uint_t faceIdOnEdge = sender->face_index(receiver->getID());
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  stencilDirection dirCellGray;
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dirCellGray = stencilDirection::CELL_GRAY_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dirCellGray = stencilDirection::CELL_GRAY_NE;
  }
  uint_t pos = 0;
  for(auto it = indexIterator(edgeIndexOnFace,
                              receiver->edge_orientation[edgeIndexOnFace],
                              CELL_GRAY, level_); it != indexIterator(); ++it){
    faceData[*it] = edgeData[BubbleEdge::edge_index(level_,pos,dirCellGray)];
    pos++;
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace BubbleFace;
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  for(auto it = indexIterator(edgeIndexOnFace,
                              sender->edge_orientation[edgeIndexOnFace],
                              CELL_BLUE, level_);
      it != indexIterator();
      ++it){
    buffer << faceData[*it];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  uint_t vPerEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t faceIdOnEdge = receiver->face_index(sender);
  stencilDirection dirCellBlue;
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dirCellBlue = stencilDirection::CELL_BLUE_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dirCellBlue = stencilDirection::CELL_BLUE_NW;
  }
  //unpack Blue Cell
  for (uint_t i = 1; i < vPerEdge - 1; ++i)
  {
    buffer >> edgeData[BubbleEdge::edge_index(level_,i,dirCellBlue)];
  }
}

template< typename ValueType >
void DGPackInfo< ValueType >::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::BubbleFace;
  ValueType *edgeData = receiver->getData(dataIDEdge_)->getPointer( level_ );
  ValueType *faceData = sender->getData(dataIDFace_)->getPointer( level_ );
  const uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  const uint_t edgeIndexOnFace = sender->edge_index(receiver->getID());
  uint_t pos = 1;
  stencilDirection dirCellBlue;
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dirCellBlue = stencilDirection::CELL_BLUE_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dirCellBlue = stencilDirection::CELL_BLUE_NW;
  }
  for(auto it = BubbleFace::indexIterator(edgeIndexOnFace,
                                          sender->edge_orientation[edgeIndexOnFace],
                                          CELL_BLUE,
                                          level_); it != indexIterator(); ++it)
  {
    edgeData[BubbleEdge::edge_index(level_,pos,dirCellBlue)] = faceData[*it];
    pos++;
  }
}


}//namespace hhg
