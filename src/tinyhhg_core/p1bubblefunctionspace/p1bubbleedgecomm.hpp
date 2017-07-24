#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubbleedgeindex.hpp"

namespace hhg {
namespace P1BubbleEdge {

using walberla::uint_t;

/*!
 * Packs data on the edge into the \ref sendBuffer.
 * Only packs data owned by edge but not halo data
 * @param edge Edge
 * @param sendBuffer Buffer to pack into
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 */
inline void packData(Edge &edge, uint_t memory_id, walberla::mpi::SendBuffer &sendBuffer, uint_t level, const Face &face) {
  auto& edge_data = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level];
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  //TODO change to index function
  for (uint_t i = 0; i < rowsize; ++i) {
    sendBuffer << edge_data[i];
  }
  if(edge.face_index(face) == 0){
    sendBuffer << edge_data[EdgeCoordsVertex::edge_index(level,0,EdgeCoordsVertex::CELL_GRAY_SE)];
    sendBuffer << edge_data[EdgeCoordsVertex::edge_index(level,rowsize-1,EdgeCoordsVertex::CELL_GRAY_SW)];
  } else {
    sendBuffer << edge_data[EdgeCoordsVertex::edge_index(level,0,EdgeCoordsVertex::CELL_GRAY_NE)];
    sendBuffer << edge_data[EdgeCoordsVertex::edge_index(level,rowsize-1,EdgeCoordsVertex::CELL_GRAY_NW)];
  }
}

/*!
 * Packs for \ref vertex data on the edge into the \ref sendBuffer.
 * @param edge Edge
 * @param sendBuffer Buffer to pack into
 * @param memory_id Memory id of the data
 * @param level Multigrid level
 * @param vertex Vertex to send to
 */
inline void packDataforVertex(Edge &edge, uint_t memory_id, walberla::mpi::SendBuffer &sendBuffer, uint_t level, const Vertex &vertex) {
  auto& edge_data = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[level];
  uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
  uint_t vertex_id = edge.vertex_index(vertex);
  if(vertex_id == 0){
    sendBuffer << edge_data[1];
  } else if(vertex_id == 1){
    sendBuffer << edge_data[levelinfo::num_microvertices_per_edge(level)-2];
  } else {
    WALBERLA_LOG_WARNING("Vertex: " << vertex << " is not contained in Edge: " << edge);
  }
}


template<size_t Level>
inline void unpackFaceData_tmpl(Edge &edge, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Face &face)
{
  auto &edge_data = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[Level];
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t pos = edge.face_index(face);
  EdgeCoordsVertex::DirVertex dir;
  EdgeCoordsVertex::DirVertex dirCellGray;
  EdgeCoordsVertex::DirVertex dirCellBlue;
  real_t dump;
  //the first edge is the south edge and the second the north edge
  if(pos == 0)
  {
    dir = EdgeCoordsVertex::VERTEX_SE;
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_SE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_SE;
  }
  else if(pos == 1)
  {
    dir = EdgeCoordsVertex::VERTEX_N;
    dirCellGray = EdgeCoordsVertex::CELL_GRAY_NE;
    dirCellBlue = EdgeCoordsVertex::CELL_BLUE_NW;
  }
  //unpack Vertex
  for (uint_t i = 0; i < rowsize -1; ++i)
  {
    recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(i,dir)];
  }
  //unpack Gray Cell
  recvBuffer >> dump;
  for (uint_t i = 1; i < rowsize - 2; ++i){
    recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(i,dirCellGray)];
  }
  recvBuffer >> dump;
  for (uint_t i = 0; i < rowsize - 2; ++i){
    recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(i,dirCellBlue)];
  }
}
SPECIALIZE(void, unpackFaceData_tmpl, unpackFaceData)

template<size_t Level>
inline void unpackVertexData_tmpl(Edge &edge, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Vertex &vertex)
{
  auto &edge_data = P1Bubble::getEdgeFunctionMemory(edge, memory_id)->data[Level];
  uint_t nbr_edges = vertex.edges.size();
  uint_t edge_id = vertex.edge_index(edge);
  uint_t vertex_id = edge.vertex_index(vertex);
  //position in edge memory
  uint_t pos;
  EdgeCoordsVertex::DirVertex dir1;
  EdgeCoordsVertex::DirVertex dir2;
  ///the vertex sends face[0] first wich is the south face by convention
  if(vertex_id == 0){
    dir1 = EdgeCoordsVertex::CELL_GRAY_SE;
    dir2 = EdgeCoordsVertex::CELL_GRAY_NE;
    pos = 0;
  } else if(vertex_id == 1) {
    pos = levelinfo::num_microvertices_per_edge(Level) - 1;
    dir1 = EdgeCoordsVertex::CELL_GRAY_SW;
    dir2 = EdgeCoordsVertex::CELL_GRAY_NW;
  } else {
    WALBERLA_LOG_WARNING("vertex " << vertex << " is not contained in edge")
  }
  recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(pos,EdgeCoordsVertex::VERTEX_C)];
  recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(pos,dir1)];
  if(edge.faces.size() == 2) {
    recvBuffer >> edge_data[EdgeCoordsVertex::index<Level>(pos, dir2)];
  }
}

SPECIALIZE(void, unpackVertexData_tmpl, unpackVertexData)

}// namespace P1BubbleEdg8e
}// namespace hhg
