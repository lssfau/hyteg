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
    buffer << edgeData[edge_index(level_,0u,CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      buffer << edgeData[edge_index(level_, 0u, CELL_GRAY_NE)];
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
    vertexData[southFaceIdOnVertex] = edgeData[edge_index(level_,0u,CELL_GRAY_SE)];
    if(sender->getNumHigherDimNeighbors() == 2)
    {
      //also send north face data
      const uint_t northFaceIdOnVertex = receiver->face_index(sender->neighborFaces()[1]);
      vertexData[northFaceIdOnVertex] = edgeData[edge_index(level_, 0u, CELL_GRAY_NE)];
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
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver->getID().getID() << " is not in Edge: " << sender);
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


void BubblePackInfo::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace hhg::BubbleFace;
  real_t *faceData = sender->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_); it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_); it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
}

void BubblePackInfo::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  uint_t vPerEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t edgeIdOnFace = receiver->face_index(sender);
  EdgeCoordsVertex::DirVertex dirCellGray;
  EdgeCoordsVertex::DirVertex dirCellBlue;
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_SE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_NE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_NW;
  }
  //unpack Gray Cell
  for (uint_t i = 0; i < vPerEdge - 1; ++i)
  {
    buffer >> edgeData[EdgeCoordsVertex::edge_index(level_,i,dirCellGray)];
  }
  //unpack Blue Cell
  for (uint_t i = 1; i < vPerEdge - 1; ++i)
  {
    buffer >> edgeData[EdgeCoordsVertex::edge_index(level_,i,dirCellBlue)];
  }
}

void BubblePackInfo::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::BubbleFace;
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  real_t *faceData = sender->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIdOnFace = receiver->face_index(sender->getID());
  EdgeCoordsVertex::DirVertex dirCellGray;
  EdgeCoordsVertex::DirVertex dirCellBlue;
  uint_t edgeIndexOnFace = sender->edge_index(receiver->getID());
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_SE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_NE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_NW;
  }
  uint_t pos = 0;
  //copy Gray Cell
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_GRAY, level_); it != indexIterator(); ++it)
  {
    edgeData[EdgeCoordsVertex::edge_index(level_,pos,dirCellGray)] = faceData[*it];
    pos++;
  }
  pos = 1;

  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], CELL_BLUE, level_); it != indexIterator(); ++it)
  {
    edgeData[EdgeCoordsVertex::edge_index(level_,pos,dirCellGray)] = faceData[*it];
    pos++;
  }
}

///@}

} //namespace hhg
