#include "P1PackInfo.hpp"
#include "p1edge.hpp"
#include "P1EdgeIndex.hpp"
#include "p1face.hpp"
#include "P1FaceIndex.hpp"

namespace hhg {

using namespace hhg::P1Edge;

/// @name Vertex to Edge
///@{

void P1PackInfo::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer){
  WALBERLA_UNUSED(receiver);
  real_t *vertexData = sender->getData(dataIDVertex_)->data[level_].get();
  buffer << vertexData[0];
}

void P1PackInfo::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer){
  using namespace EdgeCoordsVertex;
  real_t* edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  //position in edge memory
  uint_t pos;
  if(receiver->vertex_index(sender) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << sender.getID() << " is not in Edge: " << receiver)
  }
  buffer >> edgeData[edge_index(level_,pos,VERTEX_C)];
}

void P1PackInfo::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver){
  using namespace EdgeCoordsVertex;
  real_t *vertexData = sender->getData(dataIDVertex_)->data[level_].get();
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  uint_t pos;
  if(receiver->vertex_index(sender->getID()) == 0){
    pos = 0;
  } else if(receiver->vertex_index(sender->getID()) == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << sender << " is not in Edge: " << receiver)
  }
  edgeData[edge_index(level_,pos,VERTEX_C)] = vertexData[0];
}

///@}
/// @name Edge to Vertex
///@{

void P1PackInfo::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace EdgeCoordsVertex;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  uint_t vertexIdOnEdge = sender->vertex_index(receiver);
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    buffer << edgeData[edge_index(level_,1,VERTEX_C)];
  } else if(vertexIdOnEdge == 1){
    buffer << edgeData[edge_index(level_,levelinfo::num_microvertices_per_edge(level_)-2,VERTEX_C)];
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << receiver.getID() << " is not in Edge: " << sender);
  }
}

void P1PackInfo::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *vertexData = receiver->getData(dataIDVertex_)->data[level_].get();
  uint_t edgeIdOnVertex = receiver->edge_index(sender);
  buffer >> vertexData[edgeIdOnVertex + 1];
}

void P1PackInfo::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {
  using namespace EdgeCoordsVertex;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  real_t *vertexData = receiver->getData(dataIDVertex_)->data[level_].get();
  uint_t vertexIdOnEdge = sender->vertex_index(receiver->getID());
  uint_t edgeIdOnVertex = receiver->edge_index(sender->getID());
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    uint_t idx = edge_index(level_,1,VERTEX_C);
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else if(vertexIdOnEdge == 1){
    uint_t idx = edge_index(level_,levelinfo::num_microvertices_per_edge(level_)-2,VERTEX_C);
    vertexData[edgeIdOnVertex+1] = edgeData[idx];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << receiver << " is not contained in Edge: " << sender);
  }
}

///@}
/// @name Edge to Face
///@{

void P1PackInfo::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) {
  using namespace EdgeCoordsVertex;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);

  for (uint_t i = 0; i < v_perEdge; ++i) {
    buffer << edgeData[edge_index(level_,i,VERTEX_C)];
  }
}

void P1PackInfo::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  using namespace hhg::P1Face;
  real_t *faceData = receiver->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = receiver->edge_index(sender);
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_);
      it != indexIterator(); ++it){
    buffer >> faceData[*it];
  }
}

void P1PackInfo::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {
  using namespace hhg::P1Face;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  real_t *faceData = receiver->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = receiver->edge_index(sender->getID());
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_);
      it != indexIterator(); ++it)
  {
    faceData[*it] = edgeData[idx];
    idx++;
  }
}

///@}
/// @name Face to Edge
///@{

void P1PackInfo::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace hhg::P1Face;
  real_t *faceData = sender->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = sender->edge_index(receiver);
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], VERTEX_INNER, level_);
      it != indexIterator(); ++it){
    buffer << faceData[*it];
  }
}

void P1PackInfo::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  uint_t v_perEdge = levelinfo::num_microvertices_per_edge(level_);
  uint_t edgeIdOnFace = receiver->face_index(sender);
  EdgeCoordsVertex::DirVertex dir;
  //the first face is the south face and the second the north face
  if(edgeIdOnFace == 0)
  {
    dir = EdgeCoordsVertex::VERTEX_SE;
  }
  else if(edgeIdOnFace == 1)
  {
    dir = EdgeCoordsVertex::VERTEX_N;
  }
  for (uint_t i = 0; i < v_perEdge -1; ++i)
  {
    buffer >> edgeData[EdgeCoordsVertex::edge_index(level_,i,dir)];
  }
}

void P1PackInfo::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::P1Face;
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  real_t *faceData = sender->getData(dataIDFace_)->data[level_].get();
  uint_t faceIdOnEdge = receiver->face_index(sender->getID());
  EdgeCoordsVertex::DirVertex dir;
  uint_t edgeIdOnFace = sender->edge_index(receiver->getID());
  //the first face is the south face and the second the north face
  if(faceIdOnEdge == 0)
  {
    dir = EdgeCoordsVertex::VERTEX_SE;
  }
  else if(faceIdOnEdge == 1)
  {
    dir = EdgeCoordsVertex::VERTEX_N;
  }
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIdOnFace, sender->edge_orientation[edgeIdOnFace], VERTEX_INNER, level_);
      it != indexIterator(); ++it) {
    edgeData[hhg::P1Edge::EdgeCoordsVertex::edge_index(level_,idx,dir)] = faceData[*it];
    idx++;
  }
}

///@}

} //namespace hhg