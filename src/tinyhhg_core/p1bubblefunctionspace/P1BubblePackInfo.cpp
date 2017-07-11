#include "P1BubblePackInfo.hpp"
#include "p1bubbleedgeindex.hpp"
#include "p1bubblefaceindex.hpp"

namespace hhg {
namespace communication {

using namespace hhg::P1BubbleEdge;

/// @name Vertex to Edge
///@{

void P1BubblePackInfo::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer){
  uint_t nbr_neighbours = sender->edges.size();
  real_t *data = sender->getData(dataIDVertex_)->data[level_].get();
  //center vertex
  buffer << data[0];
  //cell data; each edge receives all cell data
  //TODO: this can be optimized such that only the data relevant for the receiver is send
  WALBERLA_UNUSED(receiver);
  for(uint_t i = 0; i < nbr_neighbours; ++i){
    buffer << data[nbr_neighbours + i];
  }
}

void P1BubblePackInfo::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer){
  real_t *data = receiver->getData(dataIDEdge_)->data[level_].get();
  hhg::Vertex *sendingVertex = storage_.lock()->getVertex(sender);
  uint_t edge_id = sendingVertex->edge_index(*receiver);
  uint_t vertex_id = receiver->vertex_index(*sendingVertex);
  //position in edge memory
  uint_t pos;
  hhg::P1BubbleEdge::EdgeCoordsVertex::DirVertex dir1;
  hhg::P1BubbleEdge::EdgeCoordsVertex::DirVertex dir2;
  if(vertex_id == 0){
    dir1 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NE;
    dir2 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SE;
    pos = 0;
  } else if(vertex_id == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
    dir1 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NW;
    dir2 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SW;
  } else {
    WALBERLA_LOG_WARNING("vertex " << sendingVertex << " is not contained in edge")
  }
  buffer >> data[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,pos,hhg::P1BubbleEdge::EdgeCoordsVertex::VERTEX_C)];
  //remove unneeded data
  real_t dump;
  if(edge_id > 1) {
    for (uint_t i = 0; i < edge_id - 1; ++i) {
      buffer >> dump;
    }
  }
  buffer >> data[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,pos,dir1)];
  buffer >> data[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,pos,dir2)];
}

void P1BubblePackInfo::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver){
  real_t *vertexData = sender->getData(dataIDVertex_)->data[level_].get();
  real_t *edgeData = receiver->getData(dataIDEdge_)->data[level_].get();
  uint_t vertexIdOnEdge = receiver->vertex_index(*sender);
  uint_t pos;
  hhg::P1BubbleEdge::EdgeCoordsVertex::DirVertex dir1;
  hhg::P1BubbleEdge::EdgeCoordsVertex::DirVertex dir2;
  if(vertexIdOnEdge == 0){
    dir1 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NE;
    dir2 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SE;
    pos = 0;
  } else if(vertexIdOnEdge == 1) {
    pos = levelinfo::num_microvertices_per_edge(level_) - 1;
    dir1 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_NW;
    dir2 = hhg::P1BubbleEdge::EdgeCoordsVertex::CELL_GRAY_SW;
  } else {
    WALBERLA_LOG_WARNING("vertex " << sender << " is not contained in edge")
  }
  edgeData[EdgeCoordsVertex::edge_index(level_,pos,EdgeCoordsVertex::VERTEX_C)] = vertexData[0];
  edgeData[EdgeCoordsVertex::edge_index(level_,pos,dir1)] = vertexData[1 + sender->edge_index(*receiver)];
  edgeData[EdgeCoordsVertex::edge_index(level_,pos,dir2)] = vertexData[1 + sender->edge_index(*receiver) + 1];
}

///@}
/// @name Edge to Vertex
///@{

void P1BubblePackInfo::packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  real_t *data = sender->getData(dataIDEdge_)->data[level_].get();
  uint_t vertexId = sender->vertex_index(*storage_.lock()->getVertex(receiver));
  //the last element would be the vertex itself so we have to send the next one
  if(vertexId == 0){
    buffer << data[1];
  } else if(vertexId == 1){
    buffer << data[levelinfo::num_microvertices_per_edge(level_)-2];
  } else {
    WALBERLA_LOG_WARNING("Vertex with ID: " << storage_.lock()->getVertex(receiver) << " is not contained in Edge: " << sender);
  }
}

void P1BubblePackInfo::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  real_t *data = receiver->getData(dataIDVertex_)->data[level_].get();
  uint_t edge_id = receiver->edge_index(*storage_.lock()->getEdge(sender));
  buffer >> data[edge_id + 1];
}

void P1BubblePackInfo::communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) {
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  real_t *vertexData = receiver->getData(dataIDVertex_)->data[level_].get();
  uint_t vertexIdOnEdge = sender->vertex_index(*receiver);
  uint_t edgeIdOnVertex = receiver->edge_index(*sender);
  //the last element would be the vertex itself so we have to send the next one
  if(vertexIdOnEdge == 0){
    vertexData[edgeIdOnVertex+1] = edgeData[1];
  } else if(vertexIdOnEdge == 1){
    vertexData[edgeIdOnVertex+1] = edgeData[levelinfo::num_microvertices_per_edge(level_)-2];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << receiver << " is not contained in Edge: " << sender);
  }
}

///@}
/// @name Edge to Face
///@{

void P1BubblePackInfo::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/, walberla::mpi::SendBuffer &buffer) {
  real_t *data = sender->getData(dataIDEdge_)->data[level_].get();
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level_);
  //TODO change to index function
  for (uint_t i = 0; i < rowsize; ++i) {
    buffer << data[i];
  }
}

void P1BubblePackInfo::unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
  using namespace hhg::P1BubbleFace;
  real_t *data = receiver->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = receiver->edge_index(*storage_.lock()->getEdge(sender));
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_); it != indexIterator(); ++it){
    buffer >> data[*it];
  }
}
void P1BubblePackInfo::communicateLocalEdgeToFace(const Edge *sender, Face *receiver) {
  using namespace hhg::P1BubbleFace;
  real_t *edgeData = sender->getData(dataIDEdge_)->data[level_].get();
  real_t *faceData = receiver->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = receiver->edge_index(*sender);
  uint_t idx = 0;
  for(auto it = indexIterator(edgeIndexOnFace, receiver->edge_orientation[edgeIndexOnFace], VERTEX, level_); it != indexIterator(); ++it){
    faceData[*it] = edgeData[idx];
    idx++;
  }
}

///@}
/// @name Face to Edge
///@{

void P1BubblePackInfo::packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) {
  using namespace hhg::P1BubbleFace;
  real_t *data = sender->getData(dataIDFace_)->data[level_].get();
  uint_t edgeIndexOnFace = sender->edge_index(*storage_.lock()->getEdge(receiver));
  for(auto it = indexIterator(edgeIndexOnFace, sender->edge_orientation[edgeIndexOnFace], VERTEX_INNER, level_); it != indexIterator(); ++it){
    buffer << data[*it];
  }
}
void P1BubblePackInfo::unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
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
void P1BubblePackInfo::communicateLocalFaceToEdge(const Face *sender, Edge *receiver) {
  using namespace hhg::P1BubbleFace;
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
    edgeData[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,idx,dir)] = faceData[*it];
    idx++;
  }
}

///@}

} //namespace communication
} //namespace hhg