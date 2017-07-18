#pragma once

#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "core/mpi/SendBuffer.h"

namespace hhg {
namespace P1BubbleVertex {

using walberla::real_t;
using walberla::uint_t;

/*!
 * Packs data which can be send to the edges.
 * @param level
 * @param vertex
 * @param memory_id
 * @param sendBuffer
 */
inline void packData(uint_t level, Vertex &vertex, uint_t memory_id, walberla::mpi::SendBuffer &sendBuffer, const Edge &edge) {
  uint_t nbrEdges = vertex.edges.size();
  auto &vertex_data = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
  //center vertex
  sendBuffer << vertex_data[0];
  Face& edgeFace0 = *edge.faces[0];
  uint_t vertexFace0idx = vertex.face_index(edgeFace0);
  sendBuffer << vertex_data[nbrEdges+vertexFace0idx];
  if(edge.faces.size() == 2) {
    Face &edgeFace1 = *edge.faces[1];
    uint_t vertexFace1idx = vertex.face_index(edgeFace1);
    sendBuffer << vertex_data[nbrEdges + vertexFace1idx];
  }
}

inline void unpackEdgeData(uint_t level, Vertex &vertex, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Edge &edge) {
  auto &vertex_data = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
  uint_t edge_id = vertex.edge_index(edge);
  recvBuffer >> vertex_data[edge_id + 1];
}

} // P1BubbleVertex
} // hhg
