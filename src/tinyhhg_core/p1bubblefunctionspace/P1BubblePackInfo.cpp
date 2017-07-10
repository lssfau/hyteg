#include "P1BubblePackInfo.hpp"
#include "p1bubbleedgeindex.hpp"

namespace hhg {
namespace communication {

void P1BubblePackInfo::packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer){
  uint_t nbr_neighbours = sender->edges.size();
  real_t *data = sender->getData(dataIDVertex_)->data[level_].get();
  //center vertex
  buffer << data[0];
  //cell data; each edge receives all cell data
  //TODO: this can be optimized such that only the data relevant for the receiver is send
  for(uint_t i = 0; i < nbr_neighbours; ++i){
    buffer << data[nbr_neighbours + i];
  }
}

void P1BubblePackInfo::unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer){
  real_t *data = receiver->getData(dataIDEdge_)->data[level_].get();
  hhg::Vertex *sendingVertex = storage_->getVertex(sender);
  uint_t nbr_edges = sendingVertex->edges.size();
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
  if(edge_id > 1) {
    for (uint_t i = 0; i < edge_id - 1; ++i) {
      auto dump = NULL;
      buffer >> dump;
    }
  }
  buffer >> data[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,pos,dir1)];
  buffer >> data[hhg::P1BubbleEdge::EdgeCoordsVertex::edge_index(level_,pos,dir2)];
}

void P1BubblePackInfo::communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver){
  //TODO
}

void P1BubblePackInfo::unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) {
//  real_t *data = receiver->getData(dataIDEdge_)->data[level_].get();
//  uint_t edge_id = storage_.gvertex.edge_index(edge);
//  recvBuffer >> data[edge_id + 1];
//}
}

void P1BubblePackInfo::packEdgeForFace(const Edge *sender, const PrimitiveID &/*receiver*/,
                                         walberla::mpi::SendBuffer &buffer) {
  real_t *data = sender->getData(dataIDEdge_)->data[level_].get();
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level_);
  //TODO change to index function
  for (uint_t i = 0; i < rowsize; ++i) {
    buffer << data[i];
  }
}



} //namespace communication
} //namespace hhg