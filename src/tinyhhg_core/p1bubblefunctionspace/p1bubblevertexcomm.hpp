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

inline void unpackEdgeData(Face &face, uint_t memory_id, walberla::mpi::RecvBuffer &recvBuffer, const Edge &edge) {
//  auto& face_data = P1Bubble::getFaceFunctionMemory(face, memory_id)->data[Level];
//  uint_t edgeIndex = face.edge_index(edge);
//  for(auto it = indexIterator(edgeIndex, face.edge_orientation[edgeIndex], VERTEX, Level); it != indexIterator(); ++it){
//    recvBuffer >> face_data[*it];
//  }
}

}
}