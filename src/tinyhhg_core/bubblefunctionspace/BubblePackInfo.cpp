#include "BubblePackInfo.hpp"
#include "BubbleEdgeIndex.hpp"
#include "BubbleFaceIndex.hpp"

namespace hhg {

using namespace hhg::BubbleEdge;

/// @name Vertex to Edge
///@{

void BubblePackInfo::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer){
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

void BubblePackInfo::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer){
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

void BubblePackInfo::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver){
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
}

///@}
/// @name Edge to Vertex
///@{

void BubblePackInfo::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace EdgeCoordsVertex;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  if(vertexIdOnEdge == 0){
    buffer << edgeData[edge_index(level_,0,CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      buffer << edgeData[edge_index(level_, 0, CELL_GRAY_NE)];
    }
  } else if(vertexIdOnEdge == 1){
    const uint_t lastPos = levelinfo::num_microvertices_per_edge(level_) -1;
    buffer << edgeData[edge_index(level_,lastPos,CELL_GRAY_SW)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      buffer << edgeData[edge_index(level_, lastPos, CELL_GRAY_NW)];
    }
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

void BubblePackInfo::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *vertexData = receiver->getData(dataIDVertex_)->data[level_].get();
  for(PrimitiveID faceId : storage_.lock()->getEdge(sender)->neighborFaces())
  {
    const uint_t faceIdOnVertex = receiver->face_index(faceId);
    buffer >> vertexData[faceIdOnVertex];
  }
}

void BubblePackInfo::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {
  using namespace EdgeCoordsVertex;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  real_t *vertexData = receiver->getData(dataIDVertex_)->data[level_].get();
  const uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());

  if(vertexIdOnEdge == 0){
    const uint_t southFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[0]);
    vertexData[southFaceIdOnVertex] = edgeData[edge_index(level_,0,CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      const uint_t northFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[1]);
      vertexData[northFaceIdOnVertex] = edgeData[edge_index(level_, 0, CELL_GRAY_NE)];
    }
  } else if(vertexIdOnEdge == 1){
    const uint_t lastPos = levelinfo::num_microvertices_per_edge(level_) -1;
    const uint_t southFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[0]);
    vertexData[southFaceIdOnVertex] = edgeData[edge_index(level_,lastPos,CELL_GRAY_SW)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      const uint_t northFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[1]);
      vertexData[northFaceIdOnVertex] = edgeData[edge_index(level_, lastPos, CELL_GRAY_NW)];
    }
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }

}


///@}
/// @name Edge to Face
///@{

void BubblePackInfo::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(buffer);
}

void BubblePackInfo::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
  WALBERLA_UNUSED(buffer);
}

void BubblePackInfo::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {
  //There is nothing to do here
  WALBERLA_UNUSED(sender);
  WALBERLA_UNUSED(receiver);
}

///@}
/// @name Face to Edge
///@{

///// DONE TILL HERE 11.08.2017 14:18 ITERATOR NEEDED /////

void BubblePackInfo::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace hhg::BubbleFace;
  real_t *data = sender->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = sender->edge_index(*storage_.lock()->getEdge(receiver));
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], VERTEX_INNER, level_); it != indexIterator(); ++it){
    buffer << data[*it];
  }
}

void BubblePackInfo::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *data = receiver->getData(dataIDEdge_)->data[level_].get();
  uint_t rowSize = levelinfo::num_microvertices_per_edge(level_);
  uint_t pos = receiver->face_index(*storage_.lock()->getFace(sender));
  EdgeCoordsVertex::DirVertex dir;
  //the first face is the south face and the second the north face
  if(pos == 0)
  {
    dir = EdgeCoordsVertex::VERTEX_SE;
  }
  else if(pos == 1)
  {
    dir = EdgeCoordsVertex::VERTEX_N;
  }
  for (uint_t i = 0; i < rowSize -1; ++i)
  {
    buffer >> data[EdgeCoordsVertex::edge_index(level_,i,dir)];
  }
}

void BubblePackInfo::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::BubbleFace;
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  real_t *faceData = sender->getData(dataIDFace_)->data[level_].get();
  uint_t facePosOnEdge = receiver->face_index(*sender);
  EdgeCoordsVertex::DirVertex dir;
  uint_t edgeIndexOnFace = sender->edge_index(*receiver);
  //the first face is the south face and the second the north face
  if(facePosOnEdge == 0)
  {
    dir = EdgeCoordsVertex::VERTEX_SE;
  }
  else if(facePosOnEdge == 1)
  {
    dir = EdgeCoordsVertex::VERTEX_N;
  }
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], VERTEX_INNER, level_); it != indexIterator(); ++it) {
    edgeData[hhg::BubbleEdge::EdgeCoordsVertex::edge_index(level_,idx,dir)] = faceData[*it];
    idx++;
  }
}

///@}

} //namespace communication
} //namespace hhg
