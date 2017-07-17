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
inline void packData(uint_t level, Vertex &vertex, uint_t memory_id, walberla::mpi::SendBuffer &sendBuffer) {
  uint_t nbr_neighbours = vertex.edges.size();
  auto &vertex_data = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
  //center vertex
  sendBuffer << vertex_data[0];
  //cell data; each edge receives all cell data
  for(uint_t i = 0; i < vertex.edges.size(); ++i){
    sendBuffer << vertex_data[nbr_neighbours + i];
  }
}

inline void unpackEdgeData(uint_t level, Vertex &vertex, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Edge &edge) {
  auto &vertex_data = P1Bubble::getVertexFunctionMemory(vertex, memory_id)->data[level];
  uint_t edge_id = vertex.edge_index(edge);
  recvBuffer >> vertex_data[edge_id + 1];
}

} // P1BubbleVertex
} // hhg
